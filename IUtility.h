#ifndef __ITENSORUTILITY_H_CMC__
#define __ITENSORUTILITY_H_CMC__
#include <vector>
#include <unordered_map>
#include "itensor/all.h"
#include "ReadInput.h"
#include "GeneralUtility.h"
using namespace std;
using namespace itensor;

// Use iprint
#ifndef iprint
#define iprint(name) myprinter(#name, (name))
#endif
void myprinter (string name, const ITensor& T)
{
    cout << name << endl;
    if (T)
        cout << "is real = " << isReal(T) << endl;
    cout << T << endl;
}

namespace iut
{
template <typename NumType, typename... Args>
inline NumType eltT (const ITensor& T, Args... args)
{
    if constexpr (is_same_v <NumType, Real>)
        return elt (T, args...);
    else
        return eltC (T, args...);
}

inline Real toReal (const ITensor& T)
{
    if (!isReal(T))
    {
        auto val = eltC(T);
        if (val.imag() < 1e-15)
        {
            return val.real();
        }
        else
        {
            cout << "Failed: " << __FUNCTION__ << ": " << val << endl;
            throw;
        }
    }
    else
    {
        return elt(T);
    }
}

ITensor copy_diagITensor (const ITensor& T, const IndexSet& inds)
{
    assert (order(T) == order(inds));
    assert (order(T) == 2);
    assert (T.inds()(1).dim() == T.inds()(2).dim());

    if (isReal (T))
    {
        vector<Real> elts;
        auto get_elt = [&elts] (Real r)
        {
            elts.push_back (r);
        };
        T.visit (get_elt);
        return diagITensor (elts, inds);
    }
    else
    {
        vector<Cplx> elts;
        auto get_elt = [&elts] (Cplx r)
        {
            elts.push_back (r);
        };
        T.visit (get_elt);
        return diagITensor (elts, inds);
    }
}

ITensor hard_copy (const ITensor& T)
{
    auto re = ITensor (T.inds());
    for(int i1 = 1; i1 <= T.inds()(1).dim(); i1++)
    {
        if (order(T) == 1)
        {
            if (isReal (T))
                re.set (i1, elt(T,i1));
            else
                re.set (i1, eltC(T,i1));
        }
        else
        {
            for(int i2 = 1; i2 <= T.inds()(2).dim(); i2++)
            {
                if (order(T) == 2)
                {
                    if (isReal (T))
                        re.set (i1, i2, elt(T,i1,i2));
                    else
                        re.set (i1, i2, eltC(T,i1,i2));
                }
                else
                {
                    for(int i3 = 1; i3 <= T.inds()(3).dim(); i3++)
                    {
                        if (order(T) == 3)
                        {
                            if (isReal (T))
                                re.set (i1, i2, i3, elt(T,i1,i2,i3));
                            else
                                re.set (i1, i2, i3, eltC(T,i1,i2,i3));
                        }
                        else
                        {
                            cout << "Order not implemented: " << order(T) << endl;
                            throw;
                        }
                    }
                }
            }
        }
    }
    return re;
}

template <typename ValType=Real>
inline ITensor Identity (const Index& ii, ValType val=1.)
{
    ITensor id (dag(ii), prime(ii));
    if constexpr (is_same_v <ValType, Cplx>)
    {
        if (abs(val.imag()) < 1e-15)
        {
            for(int i = 1; i <= ii.dim(); i++)
                id.set (i,i,val.real());
        }
        else
        {
            for(int i = 1; i <= ii.dim(); i++)
                id.set (i,i,val);
        }
    }
    else
    {
        for(int i = 1; i <= ii.dim(); i++)
            id.set (i,i,val);
    }
    return id;
}

template <typename ValType=Real>
inline ITensor Identity (const Index& i1, const Index& i2, ValType val=1.)
{
    ITensor id (i1, i2);
    for(int i = 1; i <= i1.dim(); i++)
        id.set (i,i,val);
    return id;
}

template <typename ValType=Real>
inline ITensor Identity (const IndexSet& iis, ValType val=1.)
{
    return Identity (iis(1), iis(2), val);
}

inline bool has_qn (const Index& ii)
{
    return !(nblock(ii) == 0);
}

// If the index <ii> has quantum number Nf (number of fermions) or Pf (parity),
// return the parities of each position in <ii>.
// Otherwise, return a vector of all elements 1
vector<bool> get_fermion_parity (const Index& ii)
{
    vector<bool> ps (1);
    if (nblock(ii) == 0)
    {
//        ps.resize (ii.dim()+1, 1);
//        return ps;
        cout << "Index has no quantum number" << endl;
        throw;
    }

    for(int i = 1; i <= nblock(ii); i++)
    // For each QN block
    {
        Real p = 1;
        if (qn(ii,i).hasName("Nf"))
            p = (qn(ii,i).val("Nf") % 2);
        else if (qn(ii,i).hasName("Pf"))
            p = (qn(ii,i).val("Pf") % 2);
        else
        {
            cout << "No Nf or Pf quantum number" << endl;
            throw;
        }

        for(int j = 1; j <= blocksize (ii, i); j++)
        // For each element in the block
        {
            ps.push_back (p);
        }
    }
    return ps;
}

// For fermionic tensors
inline ITensor parity_sign_tensor (const Index& ii)
{
    Index iip = prime(dag(ii));
    auto pfs = get_fermion_parity (ii);
    auto re = ITensor (ii, iip);
    for(int i = 1; i <= ii.dim(); i++)
    {
        int sign = (pfs.at(i) ? -1 : 1);
        re.set (i, i, sign);
    }
    return re;
/*
    int block_ipre = 0;
    for(int i = 1; i <= nblock(ii); i++)
    // For each QN block
    {
        Real a = 1.;
        if (qn(ii,i).hasName("Nf"))
            a = (qn(ii,i).val("Nf") % 2 == 1 ? -1. : 1.);
        else if (qn(ii,i).hasName("Pf"))
            a = (qn(ii,i).val("Pf") % 2 == 1 ? -1. : 1.);

        int bsize = blocksize (ii, i);
        for(int j = 1; j <= bsize; j++)
        // For each element in the block
        {
            int k = j + block_ipre;
            s.set (ii=k, iip=k, a);
        }
        block_ipre += bsize;
    }
    return dag(s);*/
}

ITensor SwapGate (const Index& i1_, const Index& i2_)
{
  Index i1 = dag(i1_),
        i2 = dag(i2_);
  Index i1_pr = prime(i1_),
        i2_pr = prime(i2_);

  vector<bool> pf1 = get_fermion_parity (i1),
               pf2 = get_fermion_parity (i2);

  ITensor swp (i1, i2, i1_pr, i2_pr);
  for(int i = 1; i <= i1.dim(); ++i)
      for(int j = 1; j <= i2.dim(); ++j)
      {
          if (pf1.at(i) && pf2.at(j))
              swp.set (i1=i, i2=j, i1_pr=j, i2_pr=i, -1.);
          else
              swp.set (i1=i, i2=j, i1_pr=j, i2_pr=i, 1.);
      }
  return swp;
}

Sweeps Read_sweeps (const string& fname, string key="sweeps", int nlines=std::numeric_limits<int>::max())
{
    vector<int> m, niter;
    vector<Real> cutoff, noise;

    ifstream ifs (fname);
    vector<string> lines = read_bracket (ifs, key, 0);

    auto keys = split_str<string> (lines.at(0));
    unordered_map <string, int> ii;
    for(int i = 0; i < keys.size(); i++)
        ii[keys.at(i)] = i;

    int nsweeps = 0;
    for(size_t i = 1; i < lines.size(); i++)
    {
        auto tmp = split_str<Real> (lines.at(i));
        int n = tmp.at(ii.at("nsweep"));
        nsweeps += n;
    }

    Sweeps sweeps (nsweeps);
    int isw = 1;
    if (nlines >= lines.size())
        nlines = lines.size() - 1;
    for(size_t i = 1; i <= nlines; i++)
    {
        auto tmp = split_str<Real> (lines.at(i));
        int nsweep = tmp.at(ii.at("nsweep"));
        for(int j = 0; j < nsweep; j++)
        {
            if (ii.count("minm") != 0)
                sweeps.setmindim (isw, tmp.at(ii.at("minm")));
            sweeps.setmaxdim (isw, tmp.at(ii.at("maxm")));
            sweeps.setcutoff (isw, tmp.at(ii.at("cutoff")));
            sweeps.setniter  (isw, tmp.at(ii.at("niter")));
            sweeps.setnoise  (isw, tmp.at(ii.at("noise")));
            isw++;
        }
    }
    return sweeps;
}

bool check_ortho (const ITensor& T, const Index& l, int pr=10, Real crit=1e-8)
// <l> is the open link
{
    auto Tdag = dag (prime (T, pr, l));
    auto IT = T * Tdag;
    assert (order(IT) == 2);
    assert (id(IT.inds()(1)) == id(IT.inds()(2)));

    auto I = Identity (IT.inds());
    auto d = norm(I - IT);
    if (d > crit)
    {
        cout << __FUNCTION__ << ": error = " << d << endl;
        return false;
    }
    return true;
}

// Contract <ten1> and <ten2> by the tags in <tags1> and <tags2>
void auto_contractEqual_by_tag (ITensor& ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    for(int i = 0; i < tags1.size(); i++)
    {
        auto tag1 = tags1.at(i);
        auto tag2 = tags2.at(i);
        auto ii1  = findIndex (ten1, tag1);
        auto ii2  = findIndex (ten2, tag2);
        if (ii1 != ii2)
            ten1 *= dag(delta(ii1,ii2));
    }
    ten1 *= ten2;
}

inline ITensor auto_contract_by_tag (ITensor ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    auto_contractEqual_by_tag (ten1, ten2, tags1, tags2);
    return ten1;
}

// Add <ten1> and <ten2> by the tags in <tags1> and <tags2>
void auto_addEqual_by_tag (ITensor& ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    for(int i = 0; i < tags1.size(); i++)
    {
        auto tag1 = tags1.at(i);
        auto tag2 = tags2.at(i);
        auto ii1  = findIndex (ten1, tag1);
        auto ii2  = findIndex (ten2, tag2);
        ten1.replaceInds ({ii1},{ii2});
    }
    ten1 += ten2;
}

// Add <ten1> and <ten2> by the tags in <tags1> and <tags2>
inline ITensor auto_add_by_tag (ITensor ten1, const ITensor& ten2, const vector<string>& tags1, const vector<string>& tags2)
{
    auto_addEqual_by_tag (ten1, ten2, tags1, tags2);
    return ten1;
}

inline int dim (const ITensor& T)
{
    int d = 1;
    for(const auto& ii : T.inds())
        d *= ii.dim();
    return d;
}

vector<ITensor> get_REs (const MPS& mps)
{
    int N = length (mps);
    vector<ITensor> Rs (N+1);
    Rs.at(N) = ITensor(1);
    for(int i = N-1; i >= 1; i--)
    {
        Rs.at(i) = Rs.at(i+1) * mps(i+1) * dag(prime(mps(i+1), "Link"));
    }
    return Rs;
}

vector<ITensor> get_LEs (const MPS& mps)
{
    int N = length (mps);
    vector<ITensor> Ls (N+1);
    Ls.at(1) = ITensor(1);
    for(int i = 2; i <= N; i++)
    {
        Ls.at(i) = Ls.at(i-1) * mps(i-1) * dag(prime(mps(i-1), "Link"));
    }
    return Ls;
}

// ==========
// Contract the transfer matrix

inline void contract_transfer (ITensor& E, const ITensor& A)
{
    E *= A;
    E *= dag(prime(A,"Link"));
}

inline void contract_transfer (ITensor& E, const ITensor& A, const ITensor& op)
{
    E *= A;
    E *= op;
    E.noPrime ("Site");
    E *= dag(prime(A,"Link"));
}

inline void contract_transfer (ITensor& E, const MPS& mps, int i)
{
    E *= mps(i);
    E *= dag(prime(mps(i),"Link"));
}

// Contract the transfer matrix
template <typename SitesT>
inline void contract_transfer (ITensor& E, const MPS& mps, int i, const SitesT& sites, string op)
{
    E *= mps(i);
    E *= sites.op(op,i);
    E.noPrime ("Site");
    E *= dag(prime(mps(i),"Link"));
}

inline ITensor merge_onsite_operators (const SiteSet& sites, int i, const vector<string>& ops)
{
    ITensor op_all (1.);
    for(int j = ops.size()-1; j >= 0; j--)
    {
        const auto& op = ops.at(j);
        if (op == "") continue;
        op_all *= prime (sites.op(op,i));
        op_all.mapPrime(1,0);
        op_all.mapPrime(2,1);
    }
    return op_all;
}

inline void apply_onsite_ops (ITensor& A, const SiteSet& sites, int i, const vector<string>& ops)
{
    assert (hasIndex (A, sites(i)));

    ITensor op_all = merge_onsite_operators (sites, i, ops);
    A *= op_all;
    A.mapPrime(1,0);
}

inline void apply_fermionic_sign (ITensor& A, const Index& iL)
{
    // Fermionic sign
    auto F = parity_sign_tensor (iL);
    A *= dag(F);
    A.noPrime();
}

void contract_transfer_matrix (ITensor& re, const SiteSet& sites, const MPS& mps, int i, const vector<string>& ops, Direction close, bool fermionic)
// The operators are applied in the reversed order in <ops>
{
    assert (close != BothDir);

    ITensor A = mps(i);
    apply_onsite_ops (A, sites, i, ops);
    if (fermionic)
    {
        if (ops.size() % 2 != 0 && i != 1)
        {
            // Fermionic sign
            auto iL = leftLinkIndex (mps, i);
            apply_fermionic_sign (A, iL);
        }
    }

    ITensor Ap = dag(mps(i));
    if (close == NoDir)
        Ap.prime ("Link");
    else if (close == Fromleft)
        Ap.prime (rightLinkIndex (mps, i));
    else if (close == Fromright)
        Ap.prime (leftLinkIndex (mps, i));

    re *= A;
    re *= Ap;
}
}
#endif
