#ifndef __SpecialMixedSiteSet_H_CMC__
#define __SpecialMixedSiteSet_H_CMC__
#include "itensor/all.h"
#include "SpecialFermion.h"
#include "SpecialBoson.h"
#include "ContainerUtility.h"
using namespace itensor;

class MixedBasis : public SiteSet
{
    int _maxOcc;

    public:

    int maxOcc () const { return _maxOcc; }

    MixedBasis() {}

    MixedBasis (int N, int iL, int iR, int iC, Args const& args=Args::global())
    {
        _maxOcc = args.getInt("MaxOcc");
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            bool in_scatter = (j >= iL and j <= iR);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite ({args,"SiteNumber=",j,"in_scatter",in_scatter}));
        }
        SiteSet::init(std::move(sites));
    }

    MixedBasis (int N, const vector<int>& scatter_sites, int iC, Args const& args=Args::global())
    {
        _maxOcc = args.getInt("MaxOcc");
        auto sites = SiteStore(N);
        for(int j = 1; j <= N; ++j)
        {
            bool in_scatter = vectool::in_vector (scatter_sites, j);
            if(j == iC) sites.set (j,SpecialBosonSite   ({args,"SiteNumber=",j}));
            else        sites.set (j,SpecialFermionSite ({args,"SiteNumber=",j,"in_scatter",in_scatter}));
        }
        SiteSet::init(std::move(sites));
    }

    MixedBasis (IndexSet const& is, Args const& args=Args::global())
    {
        int N = is.length();
        auto sites = SiteStore(N);
        for(auto j : range1(N))
        {
            auto ii = is(j);
            mycheck (hasTags(ii,"Boson") or hasTags(ii,"Fermion"), "unknown site index");
            if (hasTags(ii,"Boson"))
            {
                sites.set(j,SpecialBosonSite(ii,args));
                int d = dim(ii);
                _maxOcc = (d-1)/2;
            }
            else
                sites.set(j,SpecialFermionSite(ii));
        }
        SiteSet::init(std::move(sites));
    }
};
#endif
