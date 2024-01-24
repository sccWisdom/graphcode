#pragma once

#include <vector>

#include "util.hpp"

// ValueType - 存储堆中值的类型
// KeyType - 存储堆中键的类型
// IdxType - 索引堆中元素的索引类型，默认为vid_t
template <typename ValueType, typename KeyType, typename IdxType = vid_t>
/// 实现最小堆功能
class MinHeap
{
  private:
    IdxType n; // 堆中元素的数量
    // 存储堆中元素的向量，每个元素是一个键值对
    std::vector<std::pair<ValueType, KeyType>> heap;
    std::vector<IdxType> key2idx; // 将键映射到元素索引的向量

  public:
    MinHeap() : n(0), heap(), key2idx() {}

    /// @brief  当前索引为'cur'的元素向上移动
    /// @param cur
    /// @return 最终元素索引
    IdxType shift_up(IdxType cur)
    {
        if (cur == 0)
            return 0;
        // 找到父节点的索引
        IdxType p = (cur - 1) / 2;

        // 如果当前节点的值 < 父节点的值，则交换值和索引；
        // 然后递归调用shift_up将节点继续向上移动，保持最小堆性质
        if (heap[cur].first < heap[p].first) {
            std::swap(heap[cur], heap[p]);
            std::swap(key2idx[heap[cur].second], key2idx[heap[p].second]);
            return shift_up(p);
        }
        return cur;
    }

    /// @brief 将当前索引为 cur 的元素向下移动
    /// @param cur
    void shift_down(IdxType cur)
    {
        IdxType l = cur * 2 + 1; // 左孩子节点
        IdxType r = cur * 2 + 2; // 右孩子节点

        if (l >= n)
            return;

        IdxType m = cur;
        if (heap[l].first < heap[cur].first)
            m = l;
        if (r < n && heap[r].first < heap[m].first)
            m = r;

        if (m != cur) {
            std::swap(heap[cur], heap[m]);
            std::swap(key2idx[heap[cur].second], key2idx[heap[m].second]);
            shift_down(m);
        }
    }

    /// @brief 插入新的（键，值）
    /// @param value
    /// @param key
    void insert(ValueType value, KeyType key)
    {
        heap[n] = std::make_pair(value, key); // 插入值
        key2idx[key] = n++;                   // 更新索引
        IdxType cur = shift_up(n - 1);
        shift_down(cur);
    }

    /// @brief 检查堆中是否包含给定的键
    /// @param key
    /// @return
    bool contains(KeyType key)
    {
        return key2idx[key] < n && heap[key2idx[key]].second == key;
    }

    /// @brief 将给定键对应的值减少 d
    /// @param key
    /// @param d
    void decrease_key(KeyType key, ValueType d = 1)
    {
        if (d == 0)
            return;
        IdxType cur = key2idx[key];
        CHECK(cur < n && heap[cur].second == key) << "key not found";

        CHECK_GE(heap[cur].first, d) << "value cannot be negative";
        heap[cur].first -= d;
        shift_up(cur);
    }

    /// @brief 从堆中删除给定的键
    /// @param key
    /// @return
    bool remove(KeyType key)
    {
        IdxType cur = key2idx[key];
        if (cur >= n || heap[cur].second != key)
            return false;

        n--;
        if (n > 0) {
            heap[cur] = heap[n];
            key2idx[heap[cur].second] = cur;
            cur = shift_up(cur);
            shift_down(cur);
        }
        return true;
    }

    /// @brief 获取堆中最小的值和对应的键
    /// @param value
    /// @param key
    /// @return
    bool get_min(ValueType &value, KeyType &key)
    {
        if (n > 0) {
            value = heap[0].first;
            key = heap[0].second;
            return true;
        } else
            return false;
    }

    /// @brief 为堆预留足够的内存空间来容纳 nelements 个元素
    /// @param nelements
    void reserve(IdxType nelements)
    {
        n = 0;
        heap.resize(nelements);
        key2idx.resize(nelements);
    }

    /// @brief 清空堆中元素
    void clear() { n = 0; }
};
