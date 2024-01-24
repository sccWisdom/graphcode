#include "conversions.hpp"
#include "hdrf.hpp"
#include "util.hpp"
#include <utility>

HdrfPartitioner::HdrfPartitioner(std::string basefilename)
{
    Timer convert_timer;
    convert_timer.start();
    convert(basefilename, new Converter(basefilename));
    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time(); // 记录转换时间

    total_time.start();
    LOG(INFO) << "initializing partitioner";

    // fin = open(binedgelist_name(basefilename).c_str(), O_RDONLY,
    // (mode_t)0600); PCHECK(fin != -1) << "Error opening file for read";
    // // struct stat - 用于存储文件状态信息的数据结构
    // struct stat fileInfo = {0}; // 确保所有成员变量初始值为零
    // PCHECK(fstat(fin, &fileInfo) != -1) << "Error getting the file size";
    // PCHECK(fileInfo.st_size != 0) << "Error: file is empty";
    // LOG(INFO) << "file size: " << fileInfo.st_size;

    // // 将文件映射到内存中，为了快速访问文件内容。
    // // 0 - 表示映射的起始地址，这里设置为0，表示由系统自动选择合适的地址
    // // fileInfo.st_size - 表示映射的文件长度
    // //
    // PROT_READ-第三个参数表示内存区域的保护模式,'PROT_READ'表示内存区域只能进行读操作
    // //
    // MAP_SHARED-第四个参数表示共享映射方式,'MAP_SHARED'表示映射的内存可以被其他进程共享
    // // fin - 第五个参数是文件描述符，用于指定要映射的文件
    // // 0 -表示映射的文件偏移量，这里设置为0，表示从文件的起始位置开始映射
    // fin_map = (char *)mmap(0, fileInfo.st_size, PROT_READ, MAP_SHARED, fin,
    // 0);
    // // 成功映射到内存中，返回起始地址；映射失败，返回MAP_FAILED
    // if (fin_map == MAP_FAILED) {
    //     close(fin);
    //     PLOG(FATAL) << "error mapping the file";
    // }

    // filesize = fileInfo.st_size; // 文件总大小
    // fin_ptr = fin_map;           // 文件的起始内存映射地址
    // fin_end = fin_map + filesize; //
    // 文件的末尾内存映射地址,标记文件的结束位置

    // // 通过使用解引用操作符‘*’来获取指针位置的值
    // num_vertices = *(vid_t *)fin_ptr;
    // fin_ptr += sizeof(vid_t);
    // num_edges = *(size_t *)fin_ptr;
    // fin_ptr += sizeof(size_t);

    // LOG(INFO) << "num_vertices: " << num_vertices
    //           << ", num_edges: " << num_edges;

    std::ifstream fin(binedgelist_name(basefilename),
                      std::ios::binary | std::ios::ate);
    auto filesize = fin.tellg();
    LOG(INFO) << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);

    fin.read((char *)&num_vertices, sizeof(num_vertices));
    fin.read((char *)&num_edges, sizeof(num_edges));

    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;
    CHECK_EQ(sizeof(vid_t) + sizeof(size_t) + num_edges * sizeof(edge_t),
             filesize);

    pnum = FLAGS_p; // 分区数
    int num_partitions;
    num_partitions = pnum;
    // 生成[0,p-1]之间的随机整数
    dis.param(std::uniform_int_distribution<int>::param_type(0, p - 1));

    double lambda = 1.1;
    double balance_ratio = 1.05;
    max_partition_load = balance_ratio * num_edges / pnum; // edge load

    degrees.resize(num_vertices);
    num_batches = (filesize / ((std::size_t)memsize * 1024 * 1024)) + 1;
    num_edges_per_batch = (num_edges / num_batches) + 1;
    edge_load.resize(pnum);
    vertex_partition_matrix.assign(num_vertices, dense_bitset(pnum));
    true_vids.resize(num_vertices);

    part_degrees.assign(num_vertices, vector<vid_t>(pnum));
    balance_vertex_distribute.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
}

void HdrfPartitioner::batch_hdrf(vector<edge_t> &edges)
{
    for (auto &e : edges) {
        ++degrees[e.first];
        ++degrees[e.second];

        int max_p = find_max_score_partition_hdrf(e);
        update_vertex_partition_matrix(e, max_p);
        update_min_max_load(max_p);
        // save_edge(e.first,e.second,max_p);
        ++part_degrees[e.first][max_p];
        ++part_degrees[e.second][max_p];
        //        save_vertex(e.first,max_p);
        //        save_vertex(e.second,max_p);
    }
}

int HdrfPartitioner::find_max_score_partition_hdrf(edge_t &e)
{
    auto degree_u = degrees[e.first];
    auto degree_v = degrees[e.second];

    uint32_t sum;
    double max_score = 0;
    uint32_t max_p = 0;
    double bal, gv, gu;

    for (int p = 0; p < pnum; p++) {
        if (edge_load[p] >= max_partition_load) {
            continue;
        }

        gu = 0, gv = 0;
        sum = degree_u + degree_v;
        if (vertex_partition_matrix[e.first].get(p)) {
            gu = degree_u;
            gu /= sum;
            gu = 1 + (1 - gu);
        }
        if (vertex_partition_matrix[e.second].get(p)) {
            gv = degree_v;
            gv /= sum;
            gv = 1 + (1 - gv);
        }

        bal = max_load - edge_load[p];
        if (min_load != UINT64_MAX)
            bal /= epsilon + max_load - min_load;
        double score_p = gu + gv + lamda * bal;
        if (score_p < 0) {
            LOG(ERROR) << "ERROR: score_p < 0";
            LOG(ERROR) << "gu: " << gu;
            LOG(ERROR) << "gv: " << gv;
            LOG(ERROR) << "bal: " << bal;
            exit(-1);
        }
        if (score_p > max_score) {
            max_score = score_p;
            max_p = p;
        }
    }
    return max_p;
}

void HdrfPartitioner::update_min_max_load(int max_p)
{
    auto &load = ++edge_load[max_p];
    if (load > max_load)
        max_load = load;
}

void HdrfPartitioner::batch_node_assign_neighbors(vector<edge_t> &edges)
{
    for (auto &e : edges) {
        vid_t sp = balance_vertex_distribute[e.first],
              tp = balance_vertex_distribute[e.second];
        save_edge(e.first, e.second, sp);
        save_edge(e.second, e.first, tp);
    }
}

void HdrfPartitioner::update_vertex_partition_matrix(edge_t &e, int max_p)
{
    vertex_partition_matrix[e.first].set_bit_unsync(max_p);
    vertex_partition_matrix[e.second].set_bit_unsync(max_p);
    true_vids.set_bit_unsync(e.first);
    true_vids.set_bit_unsync(e.second);
}

void HdrfPartitioner::read_and_do(string opt_name)
{
    fin.seekg(sizeof(num_vertices) + sizeof(num_edges), std::ios::beg);
    std::vector<edge_t> edges;
    auto num_edges_left = num_edges;
    for (uint32_t i = 0; i < num_batches; i++) {
        auto edges_per_batch = num_edges_per_batch < num_edges_left
                                   ? num_edges_per_batch
                                   : num_edges_left;
        edges.resize(edges_per_batch);
        fin.read((char *)&edges[0], sizeof(edge_t) * edges_per_batch);
        if (opt_name == "hdrf") {
            batch_hdrf(edges);
        } else if (opt_name == "node_assignment") {
            batch_node_assign_neighbors(edges);
        } else {
            LOG(ERROR) << "no valid opt function";
        }
        num_edges_left -= edges_per_batch;
    }
}

void HdrfPartitioner::split()
{
    // is_mirrors容器，p个元素，都使用dense_bitset(num_vertices)进行初始化
    // is_mirrors - 容器中的每个元素都可以用于表示某个特定的位集合
    // std::vector<dense_bitset> is_mirrors(p, dense_bitset(num_vertices));
    // std::vector<size_t> counter(p, 0); // 计算每个分区中的edge的数量
    // vid_t capacity = num_edges * 1.05 / p + 1;
    // LOG(INFO) << "容量阈值：" << capacity;
    // while (fin_ptr < fin_end) {
    //     edge_t *e = (edge_t *)fin_ptr;
    //     fin_ptr += sizeof(edge_t);
    //
    //     dbh：根据顶点的度数决定放入哪个分区。对两个顶点中度数低的点进行哈希函数计算映射的分区
    // 将边划分至低度点的分区以优先复制高度点，减少低度点的复制
    //     vid_t w =
    //         degrees[e->first] <= degrees[e->second] ? e->first : e->second;
    //     int bucket = w % p; // 哈希函数
    //     counter[bucket]++;
    //     // 每个分区的位集合表示该分区中存在的顶点
    //     is_mirrors[bucket].set_bit_unsync(e->first);
    //     is_mirrors[bucket].set_bit_unsync(e->second);
    // }

    // // 取消内存映射，并关闭文件
    // if (munmap(fin_map, filesize) == -1) {
    //     close(fin);
    //     PLOG(FATAL) << "Error un-mmapping the file";
    // }
    // close(fin);

    // // 计算最大占用分区中的分数
    // size_t max_occupied = *std::max_element(counter.begin(), counter.end());
    // LOG(INFO) << "balance: " << (double)max_occupied / ((double)num_edges /
    // p); size_t total_mirrors = 0; // 计算总顶点数（包含replication）
    // // rep() - 循环宏 循环变量'i'从0到'p-1'进行迭代
    // rep (i, p)
    //     total_mirrors += is_mirrors[i].popcount();
    // LOG(INFO) << "total mirrors: " << total_mirrors;
    // LOG(INFO) << "replication factor: " << (double)total_mirrors /
    // num_vertices;

    // total_time.stop();
    // LOG(INFO) << "total partition time: " << total_time.get_time();

    read_and_do("hdrf");

    // 根据结点平衡性、随机分配的重叠度以及结点的度大小来判断
    size_t total_mirrors = 0; // 计算总顶点数
    vector<vid_t> buckets(num_partitions);
    double capacity = (double)true_vids.popcount() * 1.05 / num_partitions + 1;
    rep (i, num_vertices) {
        total_mirrors += vertex_partition_matrix[i].popcount();
        double max_score = 0.0;
        vid_t which_p;
        bool unique = false;
        if (vertex_partition_matrix[i].popcount() == 1) {
            unique = true;
        }
        repv (j, num_partitions) {
            if (vertex_partition_matrix[i].get(j)) {
                //                double
                //                score=((i%num_partitions==j)?1:0)+(part_degrees[i][j]/(degrees[i]+1))+(buckets[j]<
                //                capacity?1:0);
                double score = (part_degrees[i][j] / (degrees[i] + 1)) +
                               (buckets[j] < capacity ? 1 : 0);
                if (unique) {
                    which_p = j;
                } else if (max_score < score) {
                    max_score = score;
                    which_p = j;
                }
            }
        }
        ++buckets[which_p];
        save_vertex(i, which_p);
        balance_vertex_distribute[i] = which_p;
    }
    // node_fout.close();
    repv (j, num_partitions) {
        LOG(INFO) << "each partition node count: " << buckets[j];
    }

    read_and_do("node_assignment");
    edge_fout.close();

    LOG(INFO) << "total mirrors：" << total_mirrors;
    LOG(INFO) << "replication factor：" << (double)total_mirrors / num_vertices;
    LOG(INFO) << "是否相等：" << (num_vertices == true_vids.popcount() ? 1 : 0);

    // rep(i, num_partitions) LOG(INFO) << "edges in partition " << i << ": " <<
    // edge_load[i]; LOG(INFO) << "replication factor: " <<
    // (double)total_mirrors / true_vids.popcount();

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
}