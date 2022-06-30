#ifndef __BASISVISITOR_H_CMC__
#define __BASISVISITOR_H_CMC__

namespace basis
{

struct en
{
    en (int ki) : k (ki) {}
    int k;

    template <typename BasisType>
    auto operator() (const BasisType& b)   { return b.en(k); }
};

struct mu
{
    mu (int ki) : k (ki) {}
    int k;

    template <typename BasisType>
    auto operator() (const BasisType& b)   { return b.mu(k); }
};

struct name
{
    template <typename BasisType>
    auto operator() (const BasisType& b)   { return b.name(); }
};

struct size
{
    template <typename BasisType>
    auto operator() (const BasisType& b)   { return b.size(); }
};

struct C_op
{
    C_op (int i_, bool dag_) : i (i_) , dag (dag_) {}
    int i;
    bool dag;

    template <typename BasisType>
    auto operator() (const BasisType& b)   { return b.C_op(i, dag); }
};

}
#endif
