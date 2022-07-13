#ifndef __READWRITEFILE_H_CMC__
#define __READWRITEFILE_H_CMC__
#include "itensor/all.h"
#include <variant>
#include <tuple>
#include <map>
#include <unordered_map>
using namespace std;
using namespace itensor;

namespace iut {

// -------- General --------
template <typename T>
auto write (ostream& s, const T& t)
{
    itensor::write (s, t);
}
template <typename T>
auto read (istream& s, T& t)
{
    itensor::read (s, t);
}

// --------- smart pointer ---------
template <typename T>
auto write (ostream& s, const unique_ptr<T>& p)
{
    itensor::write (s, *p);
}
template <typename T>
auto read (istream& s, unique_ptr<T>& p)
{
    T t;
    itensor::read (s, t);
    p = make_unique (t);
}

// -------- Variant --------
template<typename T, typename... Ts>
auto write (std::ostream& s, const std::variant<T, Ts...>& v)
{
    std::visit ([&s](auto&& t) { write(s,t); }, v);
}
template<typename T, typename... Ts>
auto read (istream& s, std::variant<T, Ts...>& v)
{
    std::visit ([&s](auto&& t) { read(s,t); }, v);
}

// --------- tuple ---------
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
write (ostream& s, const std::tuple<Tp...>& t)
{}
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
write (ostream& s, const std::tuple<Tp...>& t)
{
    itensor::write (s, std::get<I>(t));
    write <I+1, Tp...> (s, t);
}
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
read (istream& s, std::tuple<Tp...>& t)
{}
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
read (istream& s, std::tuple<Tp...>& t)
{
    itensor::read (s, std::get<I>(t));
    read <I+1, Tp...> (s, t);
}

//---------- pair ------------
template <typename T1, typename T2>
void write (ostream& s, const std::pair<T1,T2>& p)
{
    itensor::write (s, p.first);
    itensor::write (s, p.second);
}
template <typename T1, typename T2>
void read (istream& s, std::pair<T1,T2>& p)
{
    T1 a1;
    T2 a2;
    itensor::read (s, a1);
    itensor::read (s, a2);
    p = std::make_pair (a1, a2);
}

// --------- vector ----------
template <typename T>
void write (ostream& s, const vector<T>& x)
{
    auto size = x.size();
    itensor::write (s,size);
    for(auto const& xi : x)
        write (s, xi);
}
template <typename T>
void read (istream& s, vector<T>& x)
{
    auto size = x.size();
    itensor::read (s,size);
    x.resize(size);
    for(auto& xi : x)
        read (s, xi);
}

// ------ unordered_map ------
template <typename T1, typename T2>
void write (ostream& s, const unordered_map<T1,T2>& x)
{
    auto size = x.size();
    itensor::write (s,size);
    for(auto const& [x1,x2] : x)
    {
        write (s, x1);
        write (s, x2);
    }
}
template <typename T1, typename T2>
void read (istream& s, unordered_map<T1,T2>& x)
{
    auto size = x.size();
    itensor::read (s,size);
    for(int i = 0; i < size; i++)
    {
        T1 x1;
        T2 x2;
        read (s, x1);
        read (s, x2);
        x.emplace (x1, move(x2));
    }
}

// ------ map ------
template <typename T1, typename T2>
void write (ostream& s, const map<T1,T2>& x)
{
    auto size = x.size();
    itensor::write (s,size);
    for(auto const& [x1,x2] : x)
    {
        write (s, x1);
        write (s, x2);
    }
}
template <typename T1, typename T2>
void read (istream& s, map<T1,T2>& x)
{
    auto size = x.size();
    itensor::read (s,size);
    for(int i = 0; i < size; i++)
    {
        T1 x1;
        T2 x2;
        read (s, x1);
        read (s, x2);
        x.emplace (x1, move(x2));
    }
}

}
#endif
