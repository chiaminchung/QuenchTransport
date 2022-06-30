#ifndef __CORR_H_CMC__
#define __CORR_H_CMC__
#include "itensor/all.h"
#include "IUtility.h"

template <typename SiteType, typename NumType>
Mat<NumType> exact_corr (const MPS& psi)
{
    int N = length (psi);
    auto sites = SiteType (siteInds(psi));
    auto corr = Mat<NumType> (N,N);
    for(int i = 1; i <= N; i++)
    for(int j = 1; j <= N; j++)
    {
        AutoMPO ampo (sites);
        ampo += 1.,"Cdag",i,"C",j;
        auto C = toMPO (ampo);
        if constexpr (is_same_v <NumType, Real>)
            corr(i-1,j-1) = inner (psi, C, psi);
        else
            corr(i-1,j-1) = innerC (psi, C, psi);
    }
    return corr;
}

template <typename SiteType, typename NumType>
Mat<NumType> exact_corr2 (const MPS& psi)
{
    int N = length (psi);
    auto sites = SiteType (siteInds(psi));
    auto corr = Mat<NumType> (N,N);
    for(int i = 1; i <= N; i++)
    for(int j = 1; j <= N; j++)
    {
        auto phi = psi;
        phi.ref(j) *= sites.op("C",j);
        phi.ref(j).noPrime("Site");
        phi.ref(i) *= sites.op("Cdag",i);
        phi.ref(i).noPrime("Site");
        if constexpr (is_same_v <NumType, Real>)
            corr(i-1,j-1) = inner (psi, phi);
        else
            corr(i-1,j-1) = innerC (psi, phi);
    }
    return corr;
}

inline ITensor multSite (ITensor A, ITensor const& B)
{
    A *= B;
    A.noPrime();
    return A;
}

inline Index make_index (int m)
{
    //return Index (QN({"Nf",0}), m, Out, "Mix");
    return Index (m, "Mix");
}

ITensor make_expand_proj (Index ii, int m_extra=1)
{
    ii.dag();
    // Make projector between the original space (m) and the enlarged space (m+m_extra)
    Index inew = make_index (dim(ii) + m_extra);
    auto P_expand = ITensor (ii, inew);
    for(int i = 1; i <= dim(ii); i++)
        P_expand.set (ii=i, inew=i, 1.);
    return P_expand;
}

template <typename SiteType, typename NumType>
class Cdag_Set
{
  public:
                 Cdag_Set () {}
                 Cdag_Set (const MPS& psi, int ell, string spin_str="");
    void         to_right (const MPS& psi, int ell, Args const& args=Args::global(), string spin_str="");
    Vec<NumType> apply_C  (const MPS& psi, int j, string spin_str) const;
    int          m        () const { return dim(_si); }

  private:
    vector<ITensor> _Us, _Proj;
    ITensor         _V;
    int             _ilast;
    Index           _si;	// The "site" (left) index of _V
    SiteType        _sites;
};

template <typename SiteType, typename NumType>
Cdag_Set<SiteType,NumType> :: Cdag_Set (const MPS& psi, int ell, string spin_str)
: _ilast (0)
, _sites (siteInds(psi))
{
    if (orthoCenter(psi) <= ell) {
        cout << "Error: Cdag_Set<SiteType,NumType> :: init: orthogonality center must be > ell" << endl;
        cout << "       " << orthoCenter(psi) << ", " << ell << endl;
        throw;
    }

    _V = multSite (_sites.op("Cdag"+spin_str, ell), multSite (_sites.op("F",ell), psi(ell)));
    _V *= prime (dag(psi(ell)), rightLinkIndex(psi,ell));

    // Make a dummy site-index to C
    _si = make_index (1);
    _V *= setElt(_si=1);
}

template <typename SiteType, typename NumType>
void Cdag_Set<SiteType,NumType> :: to_right (const MPS& psi, int ell, Args const& args, string spin_str)
{
    if (orthoCenter(psi) <= ell) {
        cout << "Error: Cdag_Set<SiteType,NumType> :: to_right: orthogonality center must be > ell" << endl;
        cout << "       " << orthoCenter(psi) << ", " << ell << endl;
        throw;
    }

    // Apply the i-th transfer matrix
    _V *= prime (dag(psi(ell)), "Link");
    _V *= multSite (_sites.op("F",ell), psi(ell));

    // Compute the element for c_{i=ell}
    ITensor c_l = multSite (_sites.op("Cdag"+spin_str,ell), multSite(_sites.op("F",ell), psi(ell)));
    c_l *= prime (dag(psi(ell)), rightLinkIndex(psi,ell));

    // ------- Add c_l to Cdag -------
    //
    ITensor P_ex = make_expand_proj (_si);
    _Proj.push_back (P_ex);

    _V *= P_ex; // expand _V

    // Add c_l into _V
    Index si = findIndex (_V, "Mix");
    int m = dim(si);
    _V += c_l * setElt(si=m);
    // -------------------------------------------

    // Truncate by SVD; keep the UV form
    _Us.emplace_back (si);
    ITensor Vtmp, D;
    svd (_V, _Us.back(), D, Vtmp, args);
    _V = D * Vtmp;
    _si = commonIndex (_V, _Us.back());

    _ilast++;
}

template <typename SiteType, typename NumType>
Vec<NumType> Cdag_Set<SiteType,NumType> :: apply_C (const MPS& psi, int j, string spin_str) const
{
    // Apply Cj
    ITensor V = _V * multSite (_sites.op("C"+spin_str,j), psi(j));
    V *= prime (dag(psi(j)), rightLinkIndex (psi,j-1));

    // Contract back _Us to get CidagCj
    Vec<NumType> cdagc (_ilast+1);
    for(int i = _ilast; i > 0; i--)
    {
        int im = i-1;
        V *= _Us.at(im);
        const Index& ii = V.inds()(1);
        cdagc(i) = iut::eltT<NumType> (V, ii=dim(ii));

        V *= dag(_Proj.at(im)); // Project to left space (the column of U)
    }

    Index si = V.inds()(1);
    cdagc(0) = iut::eltT<NumType> (V, si=1);

    return cdagc;
}

template <typename SiteType, typename NumType>
Mat<NumType> Measure_corr (MPS psi, Args const& args=Args::global(), string spin_str="")
{
    int i1 = args.getInt("corr_ibeg",1);
    int i2 = args.getInt("corr_iend",length(psi));
    int N = i2 - i1 + 1;

    auto corr = Mat<NumType> (N,N);
    Cdag_Set<SiteType,NumType> Cdag;

    auto sp = SiteType (siteInds(psi));
    int mmax = 0;
    for(int j = i1; j <= i2; j++)
    {
        psi.position(j);
        int jc = j - i1;

        // Diagonal
        ITensor tmp = multSite (sp.op("N"+spin_str,j), psi(j));
        ITensor ci = tmp * dag(psi(j));
        corr(jc,jc) = iut::eltT<NumType> (ci);

        if (j == i1) continue;
        else if (j == i1+1)
            Cdag = Cdag_Set<SiteType,NumType> (psi, i1, spin_str); // Initialize C
        else
            Cdag.to_right (psi, j-1, args, spin_str); // Extend C to the next site

        // Off-diagonal
        auto cdagc = Cdag.apply_C (psi, j, spin_str);
        for(int i = i1; i < j; i++)
        {
            int ic = i - i1;
            auto c = cdagc (ic);
            corr(ic,jc) = c;
            corr(jc,ic) = conjT (c);
        }

        if (Cdag.m() > mmax)
            mmax = Cdag.m();
    }
    cout << "max dim in Measure_corr = " << mmax << endl;
    return corr;
}
#endif
