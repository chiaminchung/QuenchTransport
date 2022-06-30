#ifndef __TOOLS_H_CMC__
#define __TOOLS_H_CMC__
MPO Make_NMPO (const MixedBasis& sites)
{
    AutoMPO ampo (sites);
    for(int i = 1; i <= length(sites); i++)
    {
        ampo += 1.0,"N",i;
    }
    return toMPO (ampo);
}

template <typename SiteType>
MPS make_initstate (const SiteType& sites, int Np)
{
    int N = length(sites);
    int Npar = Np;
    InitState init (sites);
    for(int i = 1; i <= N; i++)
    {
        string state;
        if (i % 2 == 0 and Np > 0)
        {
            state = "Occ";
            Np--;
        }
        else
            state = "Emp";
        init.set (i, state);
        cout << i << ": " << state << endl;
    }
    if (Np > 0)
    {
        for(int i = 1; i <= N; i += 2)
            if (Np > 0)
            {
                init.set (i,"Occ");
                Np--;
            }
    }
    auto psi = MPS (init);

    auto Nmpo = Make_NMPO (sites);
    auto Ntot = inner (psi,Nmpo,psi);
    if (Ntot != Npar)
    {
        cout << "particle number not match:" << Ntot << " " << Npar << endl;
        throw;
    }
    return psi;
}

template <typename SiteType>
MPS make_initstate_charge_site (const SiteType& sites, int Np, const WireSystem& sys)
{
    int charge_site = sys.to_glob ("C",1);
    int N = length(sites);
    int Npar = Np;
    InitState init (sites);
    int Np_device = 0;
    for(int i = 1; i <= N; i++)
    {
        string state;
        if (i != charge_site)
        {
            if (i % 2 == 0 and Np > 0)
            {
                state = "Occ";
                Np--;
                // count number of particle in device
                auto [seg,ind] = sys.to_loc (i);
                if (seg == "S")
                    Np_device++;
            }
            else
                state = "Emp";
            init.set (i, state);
            cout << i << ": " << state << endl;
        }
    }
    if (Np > 0)
    {
        for(int i = 1; i <= N; i += 2)
            if (Np > 0 and i != charge_site)
            {
                init.set (i,"Occ");
                Np--;
                cout << i << ": Occ" << endl;

                // count number of particle in device
                auto [seg,ind] = sys.to_loc (i);
                if (seg == "S")
                    Np_device++;
            }
    }
    // Set charge site
    init.set (charge_site, str(Np_device));
    cout << charge_site << ": " << Np_device << endl;
    auto psi = MPS (init);

    return psi;
}

template <typename T>
bool ordered (const vector<T>& ns, T resolution=1e-4)
{
    for(int i = 1; i < ns.size(); i++)
    {
        T n1 = ns.at(i-1),
          n2 = ns.at(i);
        if (abs(n1-n2) > resolution and n2 > n1)
            return false;
    }
    return true;
}

Real den (const SiteSet& sites, const MPS& psi, int i)
{
    AutoMPO ampo (sites);
    ampo += 1.,"N",i;
    auto NN = toMPO (ampo);
    return real(innerC(psi,NN,psi));
}

Real den (const SiteSet& sites, const MPS& psi, int i1, int i2)
{
    AutoMPO ampo (sites);
    for(int i = i1; i <= i2; i++)
        ampo += 1.,"N",i;
    auto NN = toMPO (ampo);
    return real(innerC(psi,NN,psi));
}

template <typename MPSType>
ITensor print_wf (const MPSType& psi)
{
    ITensor pp (1.);
    vector<Index> iis;
    for(int i = 1; i <= length(psi); i++)
    {
        pp *= psi(i);
        auto is = findIndex (psi(i), "Site,0");
        iis.push_back (is);
        if constexpr (is_same_v <MPO, MPSType>)
            iis.push_back (prime(is));
    }
    pp.permute (iis);
    PrintData(pp);
    return pp;
}

template <typename SiteType, typename NumType>
Mat<NumType> exact_corr2 (const MPS& psi, const WireSystem& sys)
{
    int N = length (psi);
    auto sites = SiteType (siteInds(psi));
    auto corr = Mat<NumType> (N,N);

    auto apply_op = [&sites] (MPS& phi, string op, int j)
    {
        phi.ref(j) *= sites.op(op,j);
        phi.ref(j).noPrime("Site");
    };

    int ic = sys.to_glob ("C",1);
    for(int i = 1; i <= N; i++)
        for(int j = i; j <= N; j++)
        {
            auto [seg1, iseg1] = sys.to_loc(i);
            auto [seg2, iseg2] = sys.to_loc(j);
            if (seg1 == "C" or seg2 == "C")
                continue;

            auto phi = psi;
            apply_op (phi, "C", j);
            apply_op (phi, "Cdag", i);
            if (seg1 != "S" and seg2 == "S")
                apply_op (phi, "A", ic);
            else if (seg1 == "S" and seg2 != "S")
                apply_op (phi, "Adag", ic);

            if constexpr (is_same_v <NumType, Real>)
                corr(i-1,j-1) = inner (psi, phi);
            else
                corr(i-1,j-1) = innerC (psi, phi);
            if (i != j)
                corr(j-1,i-1) = corr(i-1,j-1);
        }
    return corr;
}
#endif
