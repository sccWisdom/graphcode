#include <utility>

#include "conversions.hpp"
#include "dbh_partitioner.hpp"
#include "util.hpp"

DbhPartitioner::DbhPartitioner(std::string basefilename)
{
    Timer convert_timer;
    convert_timer.start();
    convert(basefilename, new Converter(basefilename));
    convert_timer.stop();
    LOG(INFO) << "convert time: " << convert_timer.get_time(); // 记录转换时间

    total_time.start();
    LOG(INFO) << "initializing partitioner";

    fin = open(binedgelist_name(basefilename).c_str(), O_RDONLY, (mode_t)0600);
    PCHECK(fin != -1) << "Error opening file for read";
    // struct stat - 用于存储文件状态信息的数据结构
    struct stat fileInfo = {0}; // 确保所有成员变量初始值为零
    PCHECK(fstat(fin, &fileInfo) != -1) << "Error getting the file size";
    PCHECK(fileInfo.st_size != 0) << "Error: file is empty";
    LOG(INFO) << "file size: " << fileInfo.st_size;

    // 将文件映射到内存中，为了快速访问文件内容。
    // 0 - 表示映射的起始地址，这里设置为0，表示由系统自动选择合适的地址
    // fileInfo.st_size - 表示映射的文件长度
    // PROT_READ-第三个参数表示内存区域的保护模式,'PROT_READ'表示内存区域只能进行读操作
    // MAP_SHARED-第四个参数表示共享映射方式,'MAP_SHARED'表示映射的内存可以被其他进程共享
    // fin - 第五个参数是文件描述符，用于指定要映射的文件
    // 0 -表示映射的文件偏移量，这里设置为0，表示从文件的起始位置开始映射
    fin_map = (char *)mmap(0, fileInfo.st_size, PROT_READ, MAP_SHARED, fin, 0);
    // 成功映射到内存中，返回起始地址；映射失败，返回MAP_FAILED
    if (fin_map == MAP_FAILED) {
        close(fin);
        PLOG(FATAL) << "error mapping the file";
    }

    filesize = fileInfo.st_size; // 文件总大小
    fin_ptr = fin_map;           // 文件的起始内存映射地址
    fin_end = fin_map + filesize; // 文件的末尾内存映射地址,标记文件的结束位置

    // 通过使用解引用操作符‘*’来获取指针位置的值
    num_vertices = *(vid_t *)fin_ptr;
    fin_ptr += sizeof(vid_t);
    num_edges = *(size_t *)fin_ptr;
    fin_ptr += sizeof(size_t);

    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;

    p = FLAGS_p; // 分区数
    // 生成[0,p-1]之间的随机整数
    dis.param(std::uniform_int_distribution<int>::param_type(0, p - 1));

    degrees.resize(num_vertices);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);
    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
}

void DbhPartitioner::split()
{
    // is_mirrors容器，p个元素，都使用dense_bitset(num_vertices)进行初始化
    // is_mirrors - 容器中的每个元素都可以用于表示某个特定的位集合
    std::vector<dense_bitset> is_mirrors(p, dense_bitset(num_vertices));
    std::vector<size_t> counter(p, 0); // 计算每个分区中的edge的数量
    vid_t capacity = num_edges * 1.05 / p + 1;
    LOG(INFO) << "容量阈值：" << capacity;
    while (fin_ptr < fin_end) {
        edge_t *e = (edge_t *)fin_ptr;
        fin_ptr += sizeof(edge_t);
        // dbh：根据顶点的度数决定放入哪个分区。对两个顶点中度数低的点进行哈希函数计算映射的分区
        // 将边划分至低度点的分区以优先复制高度点，减少低度点的复制
        vid_t w =
            degrees[e->first] <= degrees[e->second] ? e->first : e->second;
        int bucket = w % p; // 哈希函数
        counter[bucket]++;
        // 每个分区的位集合表示该分区中存在的顶点
        is_mirrors[bucket].set_bit_unsync(e->first);
        is_mirrors[bucket].set_bit_unsync(e->second);
    }

    // 取消内存映射，并关闭文件
    if (munmap(fin_map, filesize) == -1) {
        close(fin);
        PLOG(FATAL) << "Error un-mmapping the file";
    }
    close(fin);

    // 计算最大占用分区中的分数
    size_t max_occupied = *std::max_element(counter.begin(), counter.end());
    LOG(INFO) << "balance: " << (double)max_occupied / ((double)num_edges / p);
    size_t total_mirrors = 0; // 计算总顶点数（包含replication）
    // rep() - 循环宏 循环变量'i'从0到'p-1'进行迭代
    rep (i, p)
        total_mirrors += is_mirrors[i].popcount();
    LOG(INFO) << "total mirrors: " << total_mirrors;
    LOG(INFO) << "replication factor: " << (double)total_mirrors / num_vertices;

    total_time.stop();
    LOG(INFO) << "total partition time: " << total_time.get_time();
}
