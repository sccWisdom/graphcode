#pragma once

#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <parallel/algorithm>
#include <random>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "dense_bitset.hpp"
#include "edgepart.hpp"
#include "partitioner.hpp"
#include "util.hpp"

class HdrfPartitioner : public Partitioner
{
  private:
    const size_t BUFFER_SIZE = 64 * 1024 / sizeof(edge_t);
    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges;
    int pnum; // 分区数量
    int memsize;

    // batch processing globals
    uint32_t num_batches;
    uint32_t num_edges_per_batch;
    // double lamda;

    // use mmap for file input
    // int fin;
    // off_t filesize;
    // char *fin_map, *fin_ptr, *fin_end;

    std::vector<vid_t> degrees;
    uint64_t max_partition_load; // 最大分区负载

    // balance vertex partition distribute
    vector<vector<vid_t>> part_degrees; // each partition node degree
    // each node belongs to which unique partition
    vector<int> balance_vertex_distribute;

    vector<uint64_t> edge_load;
    vector<dense_bitset> vertex_partition_matrix;
    dense_bitset true_vids;
    uint64_t min_load = UINT64_MAX;
    uint64_t max_load = 0;
    const double epsilon = 1;

    std::random_device rd;
    // std:mt19937 伪随机数生成器类
    std::mt19937 gen;
    // std::uniform_int_distribution<int> 整数均匀分布,区间范围[]
    // std::unifrom_real_distribution<float/double> 浮点数均匀分布,区间范围[)
    std::uniform_int_distribution<int> dis;

  protected:
    void batch_hdrf(vector<edge_t> &edges);
    int find_max_score_partition_hdrf(edge_t &e);
    void update_vertex_partition_matrix(edge_t &e, int max_p);
    void update_min_max_load(int max_p);
    void batch_node_assign_neighbors(vector<edge_t> &edges);
    void read_and_do(string opt_name);

  public:
    DbhPartitioner(std::string basefilename);
    void split();
};
