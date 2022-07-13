#ifndef __HAMILTONIAN_H_CMC__
#define __HAMILTONIAN_H_CMC__

// C(i1,dag1) * C(i2,dag2) = \sum_k1 coef_i1,k1 C(k1,dag'1) * \sum_k2 coef_i2,k2 C(k2,dag'2)
// Return: vector of (coef, k1, dag'1, k2, dag'2)
template <typename Basis1, typename Basis2>
vector <tuple <auto,int,bool,int,bool>>
quadratic_operator_new (const Basis1& basis1, const Basis2& basis2, int i1, int i2, bool dag1, bool dag2, Real cutoff=1e-16)
{
    auto C1 = basis1.C_op (i1, dag1);     // i -> k, coef, dag
    auto C2 = basis2.C_op (i2, dag2);

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

template <typename Basis1, typename Basis2, typename NumType>
void add_CdagC (AutoMPO& ampo, const Basis1& basis1, const Basis2& basis2, int i1, int i2, NumType coef, const ToGlobDict& to_glob)
{
    if (i1 < 0) i1 += basis1.size() + 1;
    if (i2 < 0) i2 += basis2.size() + 1;
    vector <tuple <auto,int,bool,int,bool>> terms = quadratic_operator_new (basis1, basis2, i1, i2, true, false);

    // 
    string p1 = basis1.name(),
           p2 = basis2.name();
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
    int jc = to_glob.at({"C",1});
    for(auto [c12, k1, dag1, k2, dag2] : terms)  // coef, k1, dag1, k2, dag2
    {
        int j1 = to_glob.at({p1,k1});
        int j2 = to_glob.at({p2,k2});
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
template <typename Basis1, typename Basis2, typename NumType>
void add_SC (AutoMPO& ampo, const Basis1& basis1, const Basis2& basis2, int i1, int i2, NumType Delta, const ToGlobDict& to_glob)
{
    if (i1 < 0) i1 += basis1.size()+1;
    if (i2 < 0) i2 += basis2.size()+1;
    vector <tuple <auto,int,bool,int,bool>> terms = quadratic_operator_new (basis1, basis2, i1, i2, false, false);

    string p1 = basis1.name(),
           p2 = basis2.name();
    for(auto [c12, k1, dag1, k2, dag2] : terms)  // coef, k1, dag1, k2, dag2
    {
        int j1 = to_glob.at({p1,k1});
        int j2 = to_glob.at({p2,k2});
        if (j1 != j2)
        {
            auto c = Delta * c12;
            auto cc = iut::conj (c);
            ampo += -c, "C", j1, "C", j2;
            ampo += -cc, "Cdag", j2, "Cdag", j1;
        }
    }
}

template <typename BasisL, typename BasisR, typename BasisS, typename BasisC, typename SiteType, typename Para>
AutoMPO get_ampo_Kitaev_chain (const BasisL& leadL, const BasisR& leadR, const BasisS& scatterer, const BasisC& charge, const SiteType& sites, const Para& para, const ToGlobDict& to_glob)
{
    mycheck (length(sites) == to_glob.size(), "size not match");

    AutoMPO ampo (sites);

    // Diagonal terms
    string sname = scatterer.name();
    auto add_diag = [&ampo, &to_glob, &sname] (const auto& basis)
    {
        string p = basis.name();
        for(int i = 1; i <= basis.size(); i++)
        {
            int j = to_glob.at({p,i});
            auto en = basis.en(i);
            ampo += en, "N", j;
            if (p == sname)
            {
                auto mu = basis.mu(i);
                ampo += -0.5 * (en + mu), "I", i;
            }
        }
    };
    add_diag (leadL);
    add_diag (leadR);
    add_diag (scatterer);

    // Contact hopping
    add_CdagC (ampo, leadL, scatterer, -1, 1, -para.tcL, to_glob);
    add_CdagC (ampo, scatterer, leadL, 1, -1, -para.tcL, to_glob);
    add_CdagC (ampo, leadR, scatterer, 1, -1, -para.tcR, to_glob);
    add_CdagC (ampo, scatterer, leadR, -1, 1, -para.tcR, to_glob);

    // Charging energy
    string cname = charge.name();
    if (para.Ec != 0.)
    {
        int jc = to_glob.at({cname,1});
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
        int jc = to_glob.at({cname,1});
        ampo += para.EJ,"A2",jc;
        ampo += para.EJ,"A2dag",jc;
    }
    return ampo;
}
#endif
