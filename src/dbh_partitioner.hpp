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

class DbhPartitioner : public Partitioner
{
  private:
    const size_t BUFFER_SIZE = 64 * 1024 / sizeof(edge_t);
    std::string basefilename;

    vid_t num_vertices;
    size_t num_edges;
    int p;

    // use mmap for file input
    int fin;
    off_t filesize;
    char *fin_map, *fin_ptr, *fin_end;

    std::vector<vid_t> degrees;

    std::random_device rd;
    // std:mt19937 伪随机数生成器类
    std::mt19937 gen;
    // std::uniform_int_distribution<int> 整数均匀分布,区间范围[]
    // std::unifrom_real_distribution<float/double> 浮点数均匀分布,区间范围[)
    std::uniform_int_distribution<int> dis;

  public:
    DbhPartitioner(std::string basefilename);
    void split();
};
