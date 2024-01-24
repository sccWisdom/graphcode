#pragma once

#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "graph.hpp"
#include "min_heap.hpp"
#include "partitioner.hpp"
#include "util.hpp"

/* Streaming Neighbor Expansion (SNE) */
class SnePartitioner : public Partitioner
{
  private:
    const double BALANCE_RATIO = 1.05; // 平衡分区比例
    size_t BUFFER_SIZE;                // 存储边的缓冲区大小

    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges, assigned_edges;
    // bucket - 分区索引
    // p - 存储分区数目
    int p, bucket;
    // average_degree - 平均度数
    // local_average_degree - 局部平均度数
    double average_degree, local_average_degree;
    size_t max_sample_size;
    size_t capacity, local_capacity;

    // use mmap for file input
    int fin;
    off_t filesize;
    char *fin_map, *fin_ptr, *fin_end;

    std::vector<edge_t> buffer;
    std::vector<edge_t> sample_edges;
    graph_t adj_out, adj_in;
    MinHeap<vid_t, vid_t> min_heap;
    std::vector<size_t> occupied;
    std::vector<vid_t> degrees;
    std::vector<int8_t> master;
    std::vector<dense_bitset> is_cores, is_boundarys;
    std::vector<int8_t> results;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<vid_t> dis;

    edgepart_writer<vid_t, uint16_t> writer;

    /// @brief 检查给定edge是否可以分配到某个分区
    /// @param e
    /// @return 返回分区索引
    int check_edge(const edge_t *e)
    {
        // 遍历每个分区
        rep (i, bucket) {
            auto &is_boundary = is_boundarys[i];
            // 检查edge的两个端点是否都在 Boundary Set中
            if (is_boundary.get(e->first) && is_boundary.get(e->second) &&
                occupied[i] < capacity) {
                return i;
            }
        }

        rep (i, bucket) {
            auto &is_core = is_cores[i], &is_boundary = is_boundarys[i];
            // edge的两个顶点是否至少有一个属于该分区的Core Set
            if ((is_core.get(e->first) || is_core.get(e->second)) &&
                occupied[i] < capacity) {
                if (is_core.get(e->first) &&
                    degrees[e->second] > average_degree)
                    continue;
                if (is_core.get(e->second) &&
                    degrees[e->first] > average_degree)
                    continue;
                is_boundary.set_bit(e->first);
                is_boundary.set_bit(e->second);
                return i;
            }
        }

        return p;
    }

    void assign_edge(int bucket, vid_t from, vid_t to)
    {
        writer.save_edge(from, to, bucket);
        assigned_edges++;
        occupied[bucket]++;
        degrees[from]--;
        degrees[to]--;
    }

    void add_boundary(vid_t vid)
    {
        auto &is_core = is_cores[bucket], &is_boundary = is_boundarys[bucket];

        if (is_boundary.get(vid))
            return;
        is_boundary.set_bit_unsync(vid);

        if (!is_core.get(vid)) {
            min_heap.insert(adj_out[vid].size() + adj_in[vid].size(), vid);
        }

        rep (direction, 2) {
            adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
            for (size_t i = 0; i < neighbors.size();) {
                if (sample_edges[neighbors[i].v].valid()) {
                    vid_t &u = direction ? sample_edges[neighbors[i].v].second
                                         : sample_edges[neighbors[i].v].first;
                    if (is_core.get(u)) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid);
                        min_heap.decrease_key(vid);
                        sample_edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else if (is_boundary.get(u) &&
                               occupied[bucket] < capacity) {
                        assign_edge(bucket, direction ? vid : u,
                                    direction ? u : vid);
                        min_heap.decrease_key(vid);
                        min_heap.decrease_key(u);
                        sample_edges[neighbors[i].v].remove();
                        std::swap(neighbors[i], neighbors.back());
                        neighbors.pop_back();
                    } else
                        i++;
                } else {
                    std::swap(neighbors[i], neighbors.back());
                    neighbors.pop_back();
                }
            }
        }
    }

    // 用于将顶点和与之相关的边分配给当前分区
    void occupy_vertex(vid_t vid, vid_t d)
    {
        CHECK(!is_cores[bucket].get(vid)) << "add " << vid << " to core again";
        // 将 vid 添加到 is_cores[bucket] 位图中
        is_cores[bucket].set_bit_unsync(vid);

        if (d == 0)
            return;

        add_boundary(vid);

        for (auto &i : adj_out[vid])
            if (sample_edges[i.v].valid())
                add_boundary(sample_edges[i.v].second);
        adj_out[vid].clear();

        for (auto &i : adj_in[vid])
            if (sample_edges[i.v].valid())
                add_boundary(sample_edges[i.v].first);
        adj_in[vid].clear();
    }

    bool get_free_vertex(vid_t &vid)
    {
        // 生成一个随机的顶点索引
        vid = dis(gen);
        vid_t count = 0;
        while (count < num_vertices &&
               (adj_out[vid].size() + adj_in[vid].size() == 0 ||
                adj_out[vid].size() + adj_in[vid].size() >
                    2 * local_average_degree ||
                is_cores[bucket].get(vid))) {
            vid = (vid + ++count) % num_vertices;
        }
        if (count == num_vertices)
            return false;
        return true;
    }

    void read_more();
    void read_remaining();
    void clean_samples();
    void assign_master();
    size_t count_mirrors();

  public:
    SnePartitioner(std::string basefilename);
    void split();
};
