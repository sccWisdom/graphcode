#include "sne_partitioner.hpp"
#include "conversions.hpp"
#include "shuffler.hpp"

SnePartitioner::SnePartitioner(std::string basefilename)
    : basefilename(basefilename), rd(), gen(rd()), writer(basefilename)
{
    Timer shuffle_timer;
    shuffle_timer.start();
    convert(basefilename, new Shuffler(basefilename));
    shuffle_timer.stop();
    LOG(INFO) << "shuffle time: " << shuffle_timer.get_time();

    total_time.start();
    LOG(INFO) << "initializing partitioner";

    fin = open(shuffled_binedgelist_name(basefilename).c_str(), O_RDONLY,
               (mode_t)0600);
    PCHECK(fin != -1) << "Error opening file for read";
    struct stat fileInfo = {0};
    // 判断为假，则会执行断言
    PCHECK(fstat(fin, &fileInfo) != -1) << "Error getting the file size";
    PCHECK(fileInfo.st_size != 0) << "Error: file is empty";
    LOG(INFO) << "file size: " << fileInfo.st_size;

    fin_map = (char *)mmap(0, fileInfo.st_size, PROT_READ, MAP_SHARED, fin, 0);
    if (fin_map == MAP_FAILED) {
        close(fin);
        PLOG(FATAL) << "error mapping the file";
    }

    filesize = fileInfo.st_size;
    fin_ptr = fin_map;
    fin_end = fin_map + filesize;

    num_vertices = *(vid_t *)fin_ptr;
    fin_ptr += sizeof(vid_t);
    num_edges = *(size_t *)fin_ptr;
    fin_ptr += sizeof(size_t);

    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t),
             filesize);

    p = FLAGS_p;
    average_degree = (double)num_edges * 2 / num_vertices;
    assigned_edges = 0;
    LOG(INFO) << "inmem: " << FLAGS_inmem;
    if (!FLAGS_inmem) {
        LOG(INFO) << "sample_ratio: " << FLAGS_sample_ratio;
        max_sample_size = num_vertices * FLAGS_sample_ratio;
    } else
        max_sample_size = num_edges;
    local_average_degree = 2 * (double)max_sample_size / num_vertices;
    capacity = (double)num_edges * BALANCE_RATIO / p + 1;
    BUFFER_SIZE =
        std::min(64 * 1024 / sizeof(edge_t),
                 std::max((size_t)1, (size_t)(num_edges * 0.05 / p + 1)));
    LOG(INFO) << "buffer size: " << BUFFER_SIZE;
    occupied.assign(p, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    is_cores.assign(p, dense_bitset(num_vertices));
    is_boundarys.assign(p, dense_bitset(num_vertices));
    master.assign(num_vertices, -1);
    dis.param(
        std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
};

void SnePartitioner::read_more()
{
    // 样本中 edge 数量还没达到最大样本数 并且 文件指针 < 文件的末尾指针
    while (sample_edges.size() < max_sample_size && fin_ptr < fin_end) {
        // fin_buffer_end - 当前要读取的缓冲区的末尾
        edge_t *fin_buffer_end =
            std::min((edge_t *)fin_ptr + BUFFER_SIZE, (edge_t *)fin_end);
        size_t n = fin_buffer_end - (edge_t *)fin_ptr;
        results.resize(n);

#pragma omp parallel for
        for (size_t i = 0; i < n; i++)
            results[i] = check_edge((edge_t *)fin_ptr + i);

        for (size_t i = 0; i < n; i++) {
            edge_t *e = (edge_t *)fin_ptr + i; // 对缓冲区中的每条 edge 进行处理
            if (results[i] == p)
                sample_edges.push_back(*e);
            else
                assign_edge(results[i], e->first, e->second);
        }
        fin_ptr = (char *)fin_buffer_end;
    }

    // 用于构建邻接表
    adj_out.build(sample_edges);

    adj_in.build_reverse(sample_edges);
}

void SnePartitioner::read_remaining()
{
    auto &is_boundary = is_boundarys[p - 1], &is_core = is_cores[p - 1];

    // 如果该边尚未分配，则将其设置为边界顶点并将其分配到最后一个分区（p-1）中
    for (auto &e : sample_edges)
        if (e.valid()) {
            is_boundary.set_bit_unsync(e.first);
            is_boundary.set_bit_unsync(e.second);
            assign_edge(p - 1, e.first, e.second);
        }

    while (fin_ptr < fin_end) {
        edge_t *fin_buffer_end =
            std::min((edge_t *)fin_ptr + BUFFER_SIZE, (edge_t *)fin_end);
        size_t n = fin_buffer_end - (edge_t *)fin_ptr;
        results.resize(n);

#pragma omp parallel for
        for (size_t i = 0; i < n; i++)
            results[i] = check_edge((edge_t *)fin_ptr + i);

        for (size_t i = 0; i < n; i++) {
            edge_t *e = (edge_t *)fin_ptr + i;
            if (results[i] == p) {
                is_boundary.set_bit_unsync(e->first);
                is_boundary.set_bit_unsync(e->second);
                assign_edge(p - 1, e->first, e->second);
            } else
                assign_edge(results[i], e->first, e->second);
        }
        fin_ptr = (char *)fin_buffer_end;
    }

    repv (i, num_vertices) {
        if (is_boundary.get(i)) {
            is_core.set_bit_unsync(i);
            rep (j, p - 1)
                // 如果是其它分区的核心顶点，则将该顶点从该分区的core集中移除
                if (is_cores[j].get(i)) {
                    is_core.set_unsync(i, false);
                    break;
                }
        }
    }
}

// 确保样本中只保留未被分配的边
void SnePartitioner::clean_samples()
{
    for (size_t i = 0; i < sample_edges.size();) {
        // 样本边有效（样本边尚未分配）
        if (sample_edges[i].valid()) {
            int bucket = check_edge(&sample_edges[i]);
            if (bucket < p) {
                assign_edge(bucket, sample_edges[i].first,
                            sample_edges[i].second);
                // 移除当前边，并将最后一条边移到当前位置
                std::swap(sample_edges[i], sample_edges.back());
                sample_edges.pop_back();
            } else
                i++;
        } else {
            // 与sample_edges中末尾元素交换
            std::swap(sample_edges[i], sample_edges.back());
            // remove 交换后的末尾元素
            sample_edges.pop_back();
        }
    }
}

// 用于为 boundary 集中的顶点分配 master 分区
void SnePartitioner::assign_master()
{
    // 每个分区中分配为主分区的边界顶点数量
    std::vector<vid_t> count_master(p, 0);
    // 每个分区还可分配为主分区的边界顶点的配额
    std::vector<vid_t> quota(p, num_vertices);
    // 剩余未分配为主分区的边界顶点数量的变量
    long long sum = p * num_vertices;
    // 创建一个均匀分布的<double>分布器，用于生成随机数，范围为[0.0,1.0]
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    // 迭代器 - 用于迭代遍历每个分区的 is_boundarys(boundary集的顶点集合)
    std::vector<dense_bitset::iterator> pos(p);
    rep (b, p)
        pos[b] = is_boundarys[b].begin();
    vid_t count = 0;
    while (count < num_vertices) {
        long long r = distribution(gen) * sum;
        int k;
        for (k = 0; k < p; k++) {
            if (r < quota[k])
                break;
            r -= quota[k];
        }
        while (pos[k] != is_boundarys[k].end() && master[*pos[k]] != -1)
            pos[k]++;
        if (pos[k] != is_boundarys[k].end()) {
            count++;
            master[*pos[k]] = k;
            writer.save_vertex(*pos[k], k);
            count_master[k]++;
            quota[k]--;
            sum--;
        }
    }
    int max_masters =
        *std::max_element(count_master.begin(), count_master.end());
    LOG(INFO) << "master balance: "
              << (double)max_masters / ((double)num_vertices / p);
}

size_t SnePartitioner::count_mirrors()
{
    size_t result = 0;
    rep (i, p)
        result +=
            is_boundarys[i].popcount(); // 统计每个分区中 boundary 集的数量
    return result;
}

void SnePartitioner::split()
{
    LOG(INFO) << "partition `" << basefilename << "'";
    LOG(INFO) << "number of partitions: " << p;

    Timer read_timer, compute_timer;

    // 为 min_heap 和 sample_edges 预留空间
    min_heap.reserve(num_vertices);
    sample_edges.reserve(max_sample_size);
    LOG(INFO) << "partitioning...";
    for (bucket = 0; bucket < p - 1; bucket++) {
        std::cerr << bucket << ", ";
        read_timer.start();
        read_more();
        read_timer.stop();
        DLOG(INFO) << "sample size: " << adj_out.num_edges();
        compute_timer.start();
        local_capacity =
            FLAGS_inmem ? capacity : adj_out.num_edges() / (p - bucket);
        // 保证当前分区的负载均衡；动态选择度数最小的顶点并根据当前分区的剩余容量进行分配
        while (occupied[bucket] < local_capacity) {
            vid_t d, vid;
            if (!min_heap.get_min(d, vid)) {
                if (!get_free_vertex(vid)) {
                    DLOG(INFO)
                        << "partition " << bucket << " stop: no free vertices";
                    break;
                }
                // 计算顶点度数
                d = adj_out[vid].size() + adj_in[vid].size();
            } else {
                min_heap.remove(vid);
                /* CHECK_EQ(d, adj_out[vid].size() + adj_in[vid].size()); */
            }

            occupy_vertex(vid, d);
        }
        min_heap.clear();
        clean_samples();
        compute_timer.stop();
        LOG(INFO) << "finished part: " << bucket;
    }
    bucket = p - 1;
    std::cerr << bucket << std::endl;
    read_timer.start();
    read_remaining();
    read_timer.stop();
    LOG(INFO) << "finished part: " << bucket;

    // assign_master();
    LOG(INFO) << "expected edges in each partition: " << num_edges / p;
    rep (i, p)
        DLOG(INFO) << "edges in partition " << i << ": " << occupied[i];
    size_t max_occupied = *std::max_element(occupied.begin(), occupied.end());
    LOG(INFO) << "balance: " << (double)max_occupied / ((double)num_edges / p);
    size_t total_mirrors = count_mirrors();
    LOG(INFO) << "total mirrors: " << total_mirrors;
    LOG(INFO) << "replication factor: " << (double)total_mirrors / num_vertices;
    LOG(INFO) << "time used for graph input and construction: "
              << read_timer.get_time();
    LOG(INFO) << "time used for partitioning: " << compute_timer.get_time();

    LOG(INFO) << "delayed master assignment: ";
    assign_master();

    if (munmap(fin_map, filesize) == -1) {
        close(fin);
        PLOG(FATAL) << "Error un-mmapping the file";
    }
    close(fin);

    CHECK_EQ(assigned_edges, num_edges);

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
}
