#ifndef __MEACURRENT_H_CMC__
#define __MEACURRENT_H_CMC__
#include "MixedBasis.h"
#include "ContainerUtility.h"
//#include "MPOGen.h"
#include "MPSUtility.h"
using namespace vectool;

template <typename SiteType>
MPO get_current_mpo (const SiteType& sites, const WireSystem& sys, int i, Real cutoff=1e-16)
{
    auto [partL, i1] = get_loc (sys, i);
    auto [partR, i2] = get_loc (sys, i+1);
    AutoMPO ampo (sites);
    add_CdagC (ampo, sys, partL, partR, i1, i2, 1.);
    auto mpo = toMPO (ampo);
    return mpo;
}

template <typename SiteType>
MPO get_current_mpo (const SiteType& sites, const WireSystem& sys, const string& p1, int i1, const string& p2, int i2, Real cutoff=1e-16)
{
    const auto& chain1 = sys.parts().at(p1);
    const auto& chain2 = sys.parts().at(p2);
    if (i1 < 0) i1 += visit (basis::size(), chain1) + 1;
    if (i2 < 0) i2 += visit (basis::size(), chain2) + 1;
    AutoMPO ampo (sites);
    add_CdagC (ampo, sys, p1, p2, i1, i2, 1.);
    auto mpo = toMPO (ampo);
    return mpo;
}

template <typename SiteType>
MPO get_current_N_mpo (const SiteType& sites, const WireSystem& sys, int i, Real cutoff=1e-16)
{
    auto [partL, i1] = get_loc (sys, i);
    auto [partR, i2] = get_loc (sys, i+1);
    AutoMPO ampo (sites);
    add_CdagC (ampo, sys, partL, partR, i1, i2, 1.);
    auto mpo = toMPO (ampo);
    return mpo;
}

Real get_current (const MPO& JMPO, const MPS& psi)
{
    auto J = innerC (psi, JMPO, psi);
//cout << "J " << J << endl;
    return -2. * imag(J);
}
/*
template <typename SiteType>
MPO get_current_mpo2 (const SiteType& sites, const WireSystem& sys, int i)
{
    auto [p1, i1] = get_loc (sys, i);
    auto [p2, i2] = get_loc (sys, i+1);
    auto u1 = conj (sys.parts().at(p1).Ui(i1));
    auto u2 = sys.parts().at(p2).Ui(i2);

    int N = length(sites);
    auto gen = MPOGen (N, 2);
    for(int k1 = 1; k1 <= u1.size(); k1++)
        for(int k2 = 1; k2 <= u2.size(); k2++)
        {
            int j1 = sys.to_glob (p1,k1);
            int j2 = sys.to_glob (p2,k2);
            if (j1 == j2)
            {
                auto Nop = sites.op("N",j1);
                auto id = Identity (Nop.inds());
                auto CC = id - 2.*Nop;
                gen.set_onsite (j1, CC);
            }
            else if (j1 < j2)
            {
                gen.set_beg (j1, u1(k1-1), "Cdag", 3);
                gen.set_end (j2, u2(k2-1), "C", 3);
            }
            else if (j1 > j2)
            {
                gen.set_beg (j2, -u2(k2-1), "C", 4);
                gen.set_end (j1, u1(k1-1), "Cdag", 4);
            }
        }
    for(int k = 1; k <= N; k++)
    {
        gen.set_long_range (k, 3);
        gen.set_long_range (k, 4);
    }
    auto mpo = gen.makeMPO (sites);
    return mpo;
}*/
#endif
