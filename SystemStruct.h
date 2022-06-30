#ifndef __SYSTEMSTRUCT_H_CMC__
#define __SYSTEMSTRUCT_H_CMC__
#include <map>
#include <variant>
#include "GeneralUtility.h"
#include "StringUtility.h"
#include "ReadInput.h"
#include "ReadWriteFile.h"
#include "GateContainer.h"
#include "IUtility.h"
#include "SortBasis.h"
#include "ContainerUtility.h"
#include "BasisVisitor.h"
#include "OneParticleBasis.h"
#include "BdGBasis.h"
using namespace std;

// C(i1,dag1) * C(i2,dag2) = \sum_k1 coef_i1,k1 C(k1,dag'1) * \sum_k2 coef_i2,k2 C(k2,dag'2)
// Return: vector of (coef, k1, dag'1, k2, dag'2)
template <typename Basis1, typename Basis2>
vector <tuple <auto,int,bool,int,bool>>
quadratic_operator (const Basis1& basis1, const Basis2& basis2, int i1, int i2, bool dag1, bool dag2, Real cutoff=1e-16)
{
    auto C1 = visit (basis::C_op(i1, dag1), basis1),     // i -> k, coef, dag
         C2 = visit (basis::C_op(i2, dag2), basis2);

    vector<tuple <Real,int,bool,int,bool>> ops;             // coef, k1, dag1, k2, dag2
    for(auto&& [k1,c1,dag1p] : C1)
    {
        for(auto&& [k2,c2,dag2p] : C2)
        {
            auto coef = c1*c2;
            if (abs(coef) > cutoff)
                ops.emplace_back (coef,k1,dag1p,k2,dag2p);    // Cdag_ki1 C_ki2
        }
    }
    return ops;
}

using Basis = variant <OneParticleBasis, BdGBasis>;

class WireSystem
{
    public:
        WireSystem () {}

        void  add_partition (const string& name, const Basis& basis);

        int   N       () const;
        int   N_phys  () const { return has_part("C") ? N()-1 : N(); }
        int   idevL   () const { return visit (basis::size(), _parts.at("L")) + 1; }
        int   idevR   () const { if (has_part("C"))
                                    return visit (basis::size(), _parts.at("L")) + visit (basis::size(), _parts.at("S")) + visit (basis::size(), _parts.at("C")); 
                                 else
                                    return visit (basis::size(), _parts.at("L")) + visit (basis::size(), _parts.at("S")); }

        template <typename SortFunction, typename... FuncArgs>
        void sort_basis (SortFunction sort_func, FuncArgs... args);
        void sort_basis (const vector<BasisInfo>& orbs);


        const auto& to_glob    (const string& seg, int ki) const { return _to_glob.at(seg).at(ki); }
        const auto& to_loc     (int i)                     const { return _to_local.at(i); }
        const auto& parts      ()                          const { return _parts; }
        bool        has_part   (const string& p)           const { return _parts.count(p) != 0; }
        const auto& orbs       ()                          const { return _orbs; }
        void        print_orbs ()                          const;

        void write (ostream& s) const;
        void read  (istream& s);
        void write (const string& fname) const;
        void read  (const string& fname);

        template <typename T>
        vector<int> reorder_basis (vector<T>& ns, bool reverse=false);

    private:
        using GlobOrbDict = map <string, vector<int>>;

        map <string, Basis>         _parts;
        vector<BasisInfo>           _orbs;
        GlobOrbDict                 _to_glob;           // {partition, ki} -> ortical index
        vector<pair<string,int>>    _to_local;          // ortical index -> {partition, ki}

        void update_order ();
};

void WireSystem :: add_partition (const string& name, const Basis& basis)
{
    _parts.emplace (name, basis);
}

template <typename SortFunction, typename... FuncArgs>
void WireSystem :: sort_basis (SortFunction sort_func, FuncArgs... args)
{
    _orbs = sort_func (args...);
    update_order ();
}

void WireSystem :: sort_basis (const vector<BasisInfo>& orbs)
{
    mycheck (orbs.size() == this->N(), "size not match");
    _orbs = orbs;
    update_order ();
}

// Reorder _orbs based on <ns>, from small to large
// Update  _to_glob and _to_local
// Return the positions of swap gates
template <typename T>
vector<int> WireSystem :: reorder_basis (vector<T>& ns, bool reverse)
{
    mycheck (ns.size() == _orbs.size(), "size not match");
    vector<int> swap_pos;
    for(int i = 1; i < ns.size(); i++)
        for(int j = i; j >= 1; j--)
        {
            bool do_swap = (ns.at(j) < ns.at(j-1));
            if (reverse)
                do_swap = !do_swap;
            if (do_swap)
            {
                std::swap (ns[j], ns[j-1]);
                std::swap (_orbs[j], _orbs[j-1]);
                swap_pos.push_back (j);
            }
        }
    update_order ();
    return swap_pos;
}

// Update _to_glob and _to_local based on _orbs
void WireSystem :: update_order ()
{
    // --- Dictionary from {segment, k index} to orbital index ---
    // Note: _to_glob[seg,ki] = i,
    //       _to_local[i] = {seg,ki},
    //       both ki and i are 1-index
    _to_glob.clear();
    for(auto const& [name, chain] : _parts)
        _to_glob[name].resize (visit (basis::size(), chain)+1);
    _to_local.resize (_orbs.size()+1);
    for(int i = 1; i <= _orbs.size(); i++)
    {
        auto [seg, ki, en] = _orbs.at(i-1);
        _to_glob.at(seg).at(ki) = i;
        _to_local.at(i) = make_pair (seg, ki);
    }
}

template <typename NumType>
void add_CdagC (AutoMPO& ampo, const WireSystem& sys, const string& p1, const string& p2, int i1, int i2, NumType coef)
{
    const auto& chain1 = sys.parts().at(p1);
    const auto& chain2 = sys.parts().at(p2);
    if (i1 < 0) i1 += visit (basis::size(), chain1) + 1;
    if (i2 < 0) i2 += visit (basis::size(), chain2) + 1;
    vector <tuple <auto,int,bool,int,bool>> terms = quadratic_operator (chain1, chain2, i1, i2, true, false);
    //auto terms = CdagC_terms (chain1, chain2, i1, i2, coef);

    // 
    string op_charge = "";
    if ((p1 == "L" and p2 == "S") or    // Cdag_L C_S
        (p1 == "R" and p2 == "S"))      // Cdag_R C_S
    {
        op_charge = "A";
    }
    else if ((p1 == "S" and p2 == "L") or    // Cdag_S C_L
             (p1 == "S" and p2 == "R"))      // Cdag_S C_R
    {
        op_charge = "Adag";
    }
    // Hopping terms
    int jc = sys.to_glob ("C",1);
    for(auto [c12, k1, dag1, k2, dag2] : terms)  // coef, k1, dag1, k2, dag2
    {
        int j1 = sys.to_glob (p1,k1);
        int j2 = sys.to_glob (p2,k2);
        string op1 = (dag1 ? "Cdag" : "C");
        string op2 = (dag2 ? "Cdag" : "C");
        Real c = coef * c12;
        // hopping
        if (op_charge != "")
        {
            ampo += c, op1, j1, op_charge, jc, op2, j2;
        }
        else
        {
            ampo += c, op1, j1, op2, j2;
        }
    }
}

// Add -Delta C_i C_i+1 + h.c.
template <typename NumType>
void add_SC (AutoMPO& ampo, const WireSystem& sys, const string& p1, const string& p2, int i1, int i2, NumType Delta)
{
    const auto& chain1 = sys.parts().at(p1);
    const auto& chain2 = sys.parts().at(p2);
    if (i1 < 0) i1 += visit (basis::size(), chain1)+1;
    if (i2 < 0) i2 += visit (basis::size(), chain2)+1;
    vector <tuple <auto,int,bool,int,bool>> terms = quadratic_operator (chain1, chain2, i1, i2, false, false);

    for(auto [c12, k1, dag1, k2, dag2] : terms)  // coef, k1, dag1, k2, dag2
    {
        int j1 = sys.to_glob (p1,k1);
        int j2 = sys.to_glob (p2,k2);
        if (j1 != j2)
        {
            auto c = Delta * c12;
            auto cc = iut::conj (c);
            ampo += -c, "C", j1, "C", j2;
            ampo += -cc, "Cdag", j2, "Cdag", j1;
        }
    }
}

template <typename SiteType, typename Para>
AutoMPO get_ampo (const WireSystem& sys, const SiteType& sites, const Para& para)
{
    mycheck (length(sites) == sys.N(), "size not match");

    AutoMPO ampo (sites);
    // Diagonal terms
    for(auto const& [p, chain] : sys.parts())
    {
        for(int i = 1; i <= visit (basis::size(), chain); i++)
        {
            int j = sys.to_glob (p,i);
            auto en = visit (basis::en(i), chain);
            ampo += en, "N", j;
            if (p == "S")
            {
                auto mu = visit (basis::mu(i), chain);
                ampo += -0.5 * (en + mu), "I", i;
            }
        }
    }
    // Contact hopping
    for(auto const& [p1,p2,i1,i2,t] : para.hops)
    {
        add_CdagC (ampo, sys, p1, p2, i1, i2, -t);
        add_CdagC (ampo, sys, p2, p1, i2, i1, -t);
    }
    // Charging energy
    if (para.Ec != 0.)
    {
        int jc = sys.to_glob ("C",1);
        ampo += para.Ec,"NSqr",jc;
        ampo += para.Ec * para.Ng * para.Ng, "I", jc;
        ampo += -2.*para.Ec * para.Ng, "N", jc;
//        ampo +=  0.5*para.Ec,"NSqr",jc;
//        ampo += -0.5*para.Ec,"N",jc;
    }
    // Superconducting
    /*if (para.Delta != 0.)
    {
        auto const& chain = sys.parts().at("S");
        for(int i = 1; i < visit (basis::size(), chain); i++)
            add_SC (ampo, sys, "S", "S", i, i+1, para.Delta);
    }*/
    // Josephson hopping
    if (para.EJ != 0.)
    {
        int jc = sys.to_glob ("C",1);
        ampo += para.EJ,"A2",jc;
        ampo += para.EJ,"A2dag",jc;
    }
    return ampo;
}

inline void WireSystem :: print_orbs () const
{
    cout << "orbitals, segment, ki, energy" << endl;
    for(int i = 1; i <= _orbs.size(); i++)
    {
        auto [seg, ki, en] = _orbs.at(i-1);
        cout << i << " " << seg << " " << ki << " " << en << endl;
    }
}

int WireSystem :: N () const
{
    int L = 0;
    for(auto const& [p, chain] : _parts)
        L += visit (basis::size(), chain);
    return L;
}

void WireSystem :: write (ostream& s) const
{
    iutility::write(s,_parts);
    iutility::write(s,_orbs);
    iutility::write(s,_to_glob);
    iutility::write(s,_to_local);
}

void WireSystem :: read (istream& s)
{
    iutility::read(s,_parts);
    iutility::read(s,_orbs);
    iutility::read(s,_to_glob);
    iutility::read(s,_to_local);
}

// From global to local index in real space
tuple <string,int> get_loc (const WireSystem& sys, int i)
{
    string part = "L";
    int L_L = visit (basis::size(), sys.parts().at("L"));
    int L_S = visit (basis::size(), sys.parts().at("S"));
    if (i > L_L)
    {
        part = "S";
        i -= L_L;
        if (i > L_S)
        {
            part = "R";
            i -= L_S;
        }
    }
    return {part, i};
}

template <typename T>
void reorder_basis (MPS& psi, WireSystem& system, vector<T>& ns, const Args& args=Args::global())
{
    bool reverse = args.getBool("reverse",false);
    bool verbose = args.getBool("verbose",false);
    auto sites = siteInds (psi);
    auto gates = GateContainer();

    // Make swap gate
    auto i1 = sites(1);
    auto i2 = sites(2);
    ITensor swap (dag(i1), dag(i2), prime(i1), prime(i2));
    swap.set (1,1,1,1,1.);
    swap.set (1,2,2,1,1.);
    swap.set (2,1,1,2,1.);
    swap.set (2,2,2,2,-1.);

    gates.new_gate ("swap", swap, i1, i2);
    vector<int> swap_pos = system.reorder_basis (ns, reverse);
    for(int i : swap_pos)
        gates.add ("swap",i);
    gates.apply (psi, args);
}
#endif
