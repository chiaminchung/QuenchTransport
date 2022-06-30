#ifndef __SORTBASIS_H_CMC__
#define __SORTBASIS_H_CMC__
#include "BasisVisitor.h"

using BasisInfo = tuple <string, int, Real>; // segment name, orbital index, energy

// Input: arbitrary number of bases
// Combine all the basis states and sort by the energies
template <typename BasisType>
vector<BasisInfo> sort_by_energy (std::initializer_list<BasisType> chains)
{
    vector<BasisInfo> orbs;
    for(auto const& chain : chains)
    {
        for(int i = visit (basis::size(), chain); i > 0; i--)
            orbs.emplace_back (visit (basis::name(), chain), i, visit (basis::en(i), chain));
    }
    // Sort the orbitals based on the energies
    auto sort_func = [] (const BasisInfo& s1, const BasisInfo& s2)
    {
        return get<2>(s1) < get<2>(s2);
    };
    std::sort (orbs.begin(), orbs.end(), sort_func);
    return orbs;
}

// Sort all the basis states by energy; however put the states from <chainS> in zero energy of the other states
template <typename SysBasis, typename LeadBasis>
vector<BasisInfo>
sort_by_energy_S_middle
(const SysBasis& chainS, std::initializer_list<LeadBasis> other_chains)
{
    auto orb_S = sort_by_energy ({chainS});
    auto orbs = sort_by_energy (other_chains);
    auto it = orbs.begin();
    for(; it != orbs.end(); it++)
    {
        auto [seg,i,en] = *it;
        if (en > 0.) break;
    }
    orbs.insert (it, orb_S.begin(), orb_S.end());
    return orbs;
}

// chainS at the middle; chainC at the left of chainS
template <typename SysBasis, typename LeadBasis, typename ChargeBasis>
vector<BasisInfo>
sort_by_energy_S_middle_charging
(const SysBasis& chainS, const ChargeBasis& chainC, std::initializer_list<LeadBasis> other_chains)
{
    auto orb_C = sort_by_energy ({chainC});
    auto orbs = sort_by_energy_S_middle (chainS, other_chains);
    auto it = orbs.begin();
    for(; it != orbs.end(); it++)
    {
        auto [seg,i,en] = *it;
        if (seg == "S") break;
    }
    orbs.insert (it, orb_C.begin(), orb_C.end());
    return orbs;
}

// Put the charging site at zero energy
template <typename Basis>
vector<BasisInfo>
sort_by_energy_charging
(const Basis& chainC, std::initializer_list<Basis> other_chains)
{
    auto orb_C = sort_by_energy ({chainC});
    auto orbs = sort_by_energy (other_chains);
    auto it = orbs.begin();
    for(; it != orbs.end(); it++)
    {
        auto [seg,i,en] = *it;
        if (en > 0.) break;
    }
    orbs.insert (it, orb_C.begin(), orb_C.end());
    return orbs;
}
#endif
