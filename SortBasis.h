#ifndef __SORTBASIS_H_CMC__
#define __SORTBASIS_H_CMC__

using ToGlobDict = map <pair<string,int>, int>;       // {partition, ki} -> ortical index
using ToLocDict = vector<pair<string,int>>;         // ortical index -> {partition, ki}
using SortInfo = tuple <string, int, Real>;        // basis name, orbital index, energ

// Make dictionary to map between global and local indices.
// A local index is labelled by: {Basis name (name), Index in this basis (ki)}
// A Global index is just a number.
// Return:
//      to_glob[name,ki] = i,
//      to_local[i] = {name,ki},
//      (Both ki and i are 1-index)
tuple<ToGlobDict,ToLocDict> make_orb_dicts (const vector<SortInfo>& orbs)
{
    auto to_glob = ToGlobDict ();
    auto to_local = ToLocDict (orbs.size()+1);
    for(int i = 1; i <= orbs.size(); i++)
    {
        auto [name, ki, en] = orbs.at(i-1);
        to_glob[{name,ki}] = i;
        to_local.at(i) = make_pair (name, ki);
    }
    return {to_glob, to_local};
}

// ---- Recursive functions
void get_sort_info_recursive (const vector<SortInfo>& info) {} // Terminate the recursive function calls

template <typename BasisT, typename... Bases>
void get_sort_info_recursive (vector<SortInfo>& info, const BasisT& basis, const Bases&... bases)
{
    for(int i = basis.size(); i > 0; i--)
        info.emplace_back (basis.name(), i, basis.en(i));
    get_sort_info_recursive (info, bases...);
}

template <typename BasisT, typename... Bases>
vector<SortInfo> get_sort_info (const BasisT& basis, const Bases&... bases)
{
    vector<SortInfo> info;
    get_sort_info_recursive (info, basis, bases...);
    return info;
}
// ----

// Input: arbitrary number of bases
// Combine all the basis states and sort by the energies
template <typename... Bases>
vector<SortInfo> sort_by_energy (const Bases&... bases)
{
    vector<SortInfo> orbs = get_sort_info (bases...);
    // Sort the orbitals based on the energies
    auto sort_func = [] (const SortInfo& s1, const SortInfo& s2)
    {
        return get<2>(s1) < get<2>(s2);
    };
    std::sort (orbs.begin(), orbs.end(), sort_func);
    return orbs;
}

// Sort all the basis states by energy; however put the states from <chainS> in zero energy of the other states
template <typename SysBasis, typename LeadBasis>
vector<SortInfo>
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
vector<SortInfo>
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
template <typename BasisC, typename... Bases>
vector<SortInfo>
sort_by_energy_charging
(const BasisC& chainC, const Bases&... bases)
{
    auto orb_C = sort_by_energy (chainC);
    auto orbs = sort_by_energy (bases...);
    auto it = orbs.begin();
    for(; it != orbs.end(); it++)
    {
        auto [name,i,en] = *it;
        if (en > 0.) break;
    }
    orbs.insert (it, orb_C.begin(), orb_C.end());
    return orbs;
}
#endif
