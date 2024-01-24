#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

#include <stdlib.h>

#include "shuffler.hpp"
#include "util.hpp"

void Shuffler::init()
{
    num_vertices = 0;
    num_edges = 0;
    nchunks = 0;
    degrees.reserve(1 << 20);
    pool.setWorkerCount(2, threadpool11::Pool::Method::SYNC);
    // 计算：每个线程缓冲区能容纳多少数量的edge
    // FLAGS_memsize * 1024 * 1024 - 将其转换成字节单位
    // FLAGS_memsize * 1024 * 1024/pool.getWorkerCount() 平均分配内存给每个线程
    // FLAGS_memsize*1024*1024/pool.getWorkerCount()/sizeof(edge_t)每个线程的缓冲区可以容纳的edge_t对象的数量
    chunk_bufsize =
        FLAGS_memsize * 1024 * 1024 / pool.getWorkerCount() / sizeof(edge_t);
    // 预留内存大小
    chunk_buf.reserve(chunk_bufsize);
}

void Shuffler::finalize()
{
    // 向缓冲区中写入一个特殊的边(0,0)，并将flush参数设置为true，表示要刷新缓冲区
    cwrite(edge_t(0, 0), true);
    // 等待线程池中所有工作项完成
    pool.waitAll();
    LOG(INFO) << "num_vertices: " << num_vertices
              << ", num_edges: " << num_edges;
    LOG(INFO) << "number of chunks: " << nchunks;

    std::vector<std::ifstream> fin(nchunks);
    rep (i, nchunks) {
        fin[i].open(chunk_filename(i), std::ios::binary);
        CHECK(fin[i]) << "open chunk " << i << " failed";
    }
    std::vector<bool> finished(nchunks, false);
    // 用于写入重新排序后的edge列表
    std::ofstream fout(shuffled_binedgelist_name(basefilename),
                       std::ios::binary);
    int count = 0; // 用于统计已处理的分区
    fout.write((char *)&num_vertices, sizeof(num_vertices));
    fout.write((char *)&num_edges, sizeof(num_edges));
    while (true) {
        int i = rand() % nchunks;
        // 检查输入流是否已经达到文件末尾
        if (!fin[i].eof()) {
            edge_t e;
            fin[i].read((char *)&e, sizeof(edge_t));
            if (fin[i].eof()) {
                finished[i] = true;
                count++;
                if (count == nchunks)
                    break;
            } else
                fout.write((char *)&e, sizeof(edge_t));
        }
    }
    rep (i, nchunks)
        fin[i].close();
    fout.close();
    chunk_clean();

    fout.open(degree_name(basefilename), std::ios::binary);
    fout.write((char *)&degrees[0], num_vertices * sizeof(vid_t));
    fout.close();

    LOG(INFO) << "finished shuffle";
}

void Shuffler::add_edge(vid_t from, vid_t to)
{
    if (to == from) {
        LOG(WARNING) << "Tried to add self-edge " << from << "->" << to
                     << std::endl;
        return;
    }

    num_edges++;
    from = get_vid(from);
    to = get_vid(to);
    degrees[from]++;
    degrees[to]++;

    edge_t e(from, to);
    cwrite(e);
}

/// @brief  根据块索引构建文件名字符串
/// @param chunk
/// @return 拼接后完整的文件名字符串
std::string Shuffler::chunk_filename(int chunk)
{
    std::stringstream ss;
    ss << basefilename << "." << chunk << ".chunk";
    return ss.str();
}

/// @brief 删除对应的块文件
void Shuffler::chunk_clean()
{
    rep (i, nchunks)
        remove(chunk_filename(i).c_str());
}

void Shuffler::cwrite(edge_t e, bool flush)
{
    if (!flush)
        chunk_buf.push_back(e);
    // 如果flush为true或者chunk_buf的大小达到了设定的阈值
    if (flush || chunk_buf.size() >= chunk_bufsize) {
        work_t work;
        work.shuffler = this;
        work.nchunks = nchunks;
        // 当前的缓冲区内容与work中的缓冲区进行交换，将缓冲区内容传递给要提交的工作任务
        chunk_buf.swap(work.chunk_buf);
        // 将工作任务提交给线程池进行处理
        pool.postWork<void>(work);
        nchunks++; // 增加已处理的分区数
    }
}
