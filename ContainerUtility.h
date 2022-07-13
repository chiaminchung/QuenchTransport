#ifndef __CONTAINER_UTILITY_H_CMC__
#define __CONTAINER_UTILITY_H_CMC__
#include <vector>
#include <iostream>
using std::vector, std::ostream, std::tuple;

template <typename T>
ostream& operator<< (ostream& os, const vector<T>& v)
{
    for(const auto& i : v)
        os << i << " ";
    return os;
}

namespace vectool {

template <typename T>
vector<vector<T>> split_vector (const vector<T>& v, int itv)
{
    vector<vector<T>> re;
    int N = v.size() / itv;
    if (v.size() % 2 != 0)
        N++;
    re.reserve (N);
    for(int i = 0; i < v.size(); i += itv)
    {
        int i2 = i + itv;
        if (i2 > v.size())
            i2 = v.size();
        auto it1 = v.begin() + i;
        auto it2 = v.begin() + i2;
        re.emplace_back (it1, it2);
    }
    return re;
}

template <typename T>
inline void remove_element (vector<T>& v, const T& key)
{
    v.erase (std::remove (v.begin(), v.end(), key), v.end());
}

template <typename T>
inline bool in_vector (const vector<T>& v, const T& key)
{
    return (std::find(v.begin(), v.end(), key) != v.end());
}

template <typename T>
inline void extend_vector (vector<T>& v, const vector<T>& v2)
{
    v.insert (v.end(), v2.begin(), v2.end());
}

template <typename T, typename IteratorT>
inline tuple<bool,IteratorT> find_vector (const vector<T>& v, const T& key)
{
    auto it = std::find (std::begin(v), std::end(v), key);
    return {it != std::end(v), it};
}

template <typename T>
inline tuple<bool,int> find_vector_re_index (const vector<T>& v, const T& key)
{
    auto it = std::find (std::begin(v), std::end(v), key);
    return {it != std::end(v), it - std::begin(v)};
}

template <typename T>
inline vector<T> combine_vector (const vector<T>& v1, const vector<T>& v2)
{
    vector<T> re;
    re.reserve (v1.size() + v2.size());
    re = v1;
    re.insert (re.end(), v2.begin(), v2.end());
    return re;
}

template <typename T>
inline vector<T> sub_vector (const vector<T>& v, int i1, int i2)
{
    assert (i1 >= 0 && i1 < v.size() && i2 >= i1 && i2 < v.size());
    return vector<T> (v.begin()+i1, v.begin()+i2);
}

template <typename T>
inline void sort_vector (vector<T>& v)
{
    std::sort (v.begin(), v.end());
}

template <typename T>
inline bool check_no_duplicate (const vector<T>& v)
{
    for(int i = 0; i < v.size(); i++)
        for(int j = i+1; j < v.size(); j++)
            if (v.at(i) == v.at(j))
                return false;
    return true;
}

template <typename T>
inline bool has_common_ele (const vector<T>& ops1, const vector<T>& ops2)
{
    for(auto op1 : ops1)
        for(auto op2 : ops2)
            if (op1 == op2)
                return true;
    return false;
}

template <typename T>
inline auto min_value (const vector<T>& v)
{
    return min_element(v.begin(), v.end());
}

template <typename T>
inline auto max_value (const vector<T>& v)
{
    return max_element(v.begin(), v.end());
}

template <typename T>
inline void unordered_remove (vector<T>& v, const T& x)
{
    auto it = v.begin();
    while (it != v.end())
    {
        if (*it == x)
        {
            auto it_end = v.end() - 1;
            if (it != it_end)
            {
                *it = std::move(*it_end);
            }
            v.pop_back();
            continue;
        }
        ++it;
    }
}

} // end of namespace
#endif
