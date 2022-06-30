#include <iomanip>
#include "itensor/all.h"
#include "Timer.h"
Timers timer;
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "MPSUtility.h"
#include "MixedBasis.h"
#include "SystemStruct.h"
#include "ContainerUtility.h"
#include "Corr.h"
#include "TDVPObserver.h"
#include "MeaCurrent.h"
#include "tdvp.h"
#include "basisextension.h"
#include "InitState.h"
#include "Tools.h"
using namespace vectool;
using namespace itensor;
using namespace std;

struct Para
{
    vector<tuple<string,string,int,int,Real>> hops;
    Real Ec=0., Ng=0., Delta=0., EJ=0.;

    void write (ostream& s) const
    {
        iutility::write(s,hops);
        iutility::write(s,Ec);
        iutility::write(s,Ng);
        iutility::write(s,Delta);
        iutility::write(s,EJ);
    }

    void read (istream& s)
    {
        iutility::read(s,hops);
        iutility::read(s,Ec);
        iutility::read(s,Ng);
        iutility::read(s,Delta);
        iutility::read(s,EJ);
    }
};

void writeAll (const string& filename,
               const MPS& psi, const MPO& H,
               const WireSystem& system,
               const Para& para,
               const Args& args_basis,
               int step)
{
    ofstream ofs (filename);
    itensor::write (ofs, psi);
    itensor::write (ofs, H);
    itensor::write (ofs, args_basis);
    itensor::write (ofs, step);
    para.write (ofs);
    system.write (ofs);
}

void readAll (const string& filename,
              MPS& psi, MPO& H,
              WireSystem& system,
              Para& para,
              Args& args_basis,
              int& step)
{
    ifstream ifs = open_file (filename);
    itensor::read (ifs, psi);
    itensor::read (ifs, H);
    itensor::read (ifs, args_basis);
    itensor::read (ifs, step);
    para.read (ifs);
    system.read (ifs);
}

template <typename SiteType>
MPS apply (const SiteType& sites, MPS psi, Real coef, string opstr, int i)
{
    for(int j = 1; j <= i; j++)
    {
        auto op = (j == i ? coef*sites.op (opstr, j) : sites.op ("F", j));
        psi.ref(j) *= op;
        psi.ref(j).noPrime("Site");
    }
    psi.position(1);
    return psi;
}

template <typename SiteType>
MPS apply_ops (const SiteType& sites, MPS psi, const vector<tuple<int,auto,bool>>& ops)
{
    MPS re;
    for(auto [k, coef, dag] : ops)
    {
        string op = (dag ? "Cdag" : "C");
        auto tmp = apply (sites, psi, coef, op, k);
        if (!re)
            re = tmp;
        else if (norm(tmp) != 0.)
            re.plusEq (tmp);
    }
    return re;
}

int main(int argc, char* argv[])
{
    string infile = argv[1];
    InputGroup input (infile,"basic");

    auto L_lead   = input.getInt("L_lead");
    auto L_device   = input.getInt("L_device");
    auto t_lead     = input.getReal("t_lead");
    auto t_device   = input.getReal("t_device");
    auto t_contactL = input.getReal("t_contactL");
    auto t_contactR = input.getReal("t_contactR");
    auto mu_leadL   = input.getReal("mu_leadL");
    auto mu_leadR   = input.getReal("mu_leadR");
    auto mu_device  = input.getReal("mu_device");
    auto mu_biasL   = input.getReal("mu_biasL");
    auto mu_biasS   = input.getReal("mu_biasS");
    auto mu_biasR   = input.getReal("mu_biasR");
    auto Delta      = input.getReal("Delta");
    auto Ec         = input.getReal("Ec");
    auto Ng         = input.getReal("Ng");
    auto EJ         = input.getReal("EJ");
    auto damp_decay_length = input.getInt("damp_decay_length",0);
    auto maxCharge  = input.getInt("maxCharge");

    auto dt            = input.getReal("dt");
    auto time_steps    = input.getInt("time_steps");
    auto NumCenter     = input.getInt("NumCenter");
    auto Truncate      = input.getYesNo("Truncate");
    auto mixNumCenter  = input.getYesNo("mixNumCenter",false);
    auto globExpanNStr       = input.getString("globExpanN","inf");
    int globExpanN;
    if (globExpanNStr == "inf" or globExpanNStr == "Inf" or globExpanNStr == "INF")
        globExpanN = std::numeric_limits<int>::max();
    else
        globExpanN = std::stoi (globExpanNStr);
    auto globExpanItv        = input.getInt("globExpanItv",1);
    auto globExpanCutoff     = input.getReal("globExpanCutoff",1e-8);
    auto globExpanKrylovDim  = input.getInt("globExpanKrylovDim",3);
    auto globExpanHpsiCutoff = input.getReal("globExpanHpsiCutoff",1e-8);
    auto globExpanHpsiMaxDim = input.getInt("globExpanHpsiMaxDim",300);
    auto globExpanMethod     = input.getString("globExpanMethod","DensityMatrix");

    auto UseSVD        = input.getYesNo("UseSVD",true);
    auto SVDmethod     = input.getString("SVDMethod","gesdd");  // can be also "ITensor"
    auto WriteDim      = input.getInt("WriteDim");

    auto write         = input.getYesNo("write",false);
    auto write_dir     = input.getString("write_dir",".");
    auto write_file    = input.getString("write_file","");
    auto read          = input.getYesNo("read",false);
    auto read_dir      = input.getString("read_dir",".");
    auto read_file     = input.getString("read_file","");

    auto sweeps        = iut::Read_sweeps (infile, "sweeps");

    auto basis = input.getString("scatter_basis","SC");

    cout << setprecision(14) << endl;

    MPS psi;
    MPO H;
    // Define 
    auto glob_basis = WireSystem ();
    int step = 1;
    auto sites = MixedBasis();
    Para para;
    Args args_basis;
    // Initialization
    if (!read)
    {
        // Factor for exponential-decay hopping
        Real damp_fac = (damp_decay_length == 0 ? 1. : exp(-1./damp_decay_length));
        // Single-particle Hamiltonians
        cout << "H left lead" << endl;
        Basis leadL = OneParticleBasis ("L", L_lead, t_lead, mu_leadL, damp_fac, true, true);
        cout << "H right lead" << endl;
        Basis leadR = OneParticleBasis ("R", L_lead, t_lead, mu_leadR, damp_fac, false, true);
        cout << "H dev" << endl;
        //Basis system = OneParticleBasis ("S", L_device, t_device, mu_device, damp_fac, true, true);
        Basis system = BdGBasis ("S", L_device, t_device, mu_device, Delta);
        Basis charge = OneParticleBasis ("C", 1);

        // WireSystem
        glob_basis.add_partition ("L",leadL);
        glob_basis.add_partition ("R",leadR);
        glob_basis.add_partition ("S",system);
        glob_basis.add_partition ("C",charge);
        auto basis_info = vector<BasisInfo> ();
        if (basis == "SC")
        {
            basis_info = sort_by_energy_charging (glob_basis.parts().at("C"), {glob_basis.parts().at("L"), glob_basis.parts().at("R"), glob_basis.parts().at("S")});
        }
        else
        {
            cout << "Unknown basis: " << basis << endl;
            throw;
        }
        //glob_basis.sort_basis (sort_by_energy_S_middle (system.part("S"), {system.part("L"), system.part("R")}));
        //glob_basis.sort_basis (sort_by_energy_S_middle_charging (glob_basis.parts().at("S"), glob_basis.parts().at("C"),
        //                       {glob_basis.parts().at("L"), glob_basis.parts().at("R")}));
        glob_basis.sort_basis (basis_info);
        glob_basis.print_orbs();
        cout << "device site = " << glob_basis.idevL() << " " << glob_basis.idevR() << endl;

        // SiteSet
        int N = glob_basis.N();
        int charge_site = glob_basis.to_glob ("C",1);
        // Find sites in S
        vector<int> scatter_sites;
        for(int i = 0; i < glob_basis.orbs().size(); i++)
        {
            if (get<0>(glob_basis.orbs().at(i)) == "S")
                scatter_sites.push_back (i+1);
        }
        // Make SiteSet
        string systype = "Normal";
        if (Delta != 0.)
        {
            if (EJ != 0.)
                systype = "SC_Josephson_scatter";
            else
                systype = "SC_scatter";
        }
        args_basis = {"MaxOcc",maxCharge,"SystemType",systype};
        sites = MixedBasis (N, scatter_sites, charge_site, args_basis);
        cout << "charge site = " << charge_site << endl;

        // Make Hamiltonian MPO
        para.Ec = Ec;   para.Ng = Ng;   para.Delta = Delta;  para.EJ = EJ;
        para.hops.clear();
        para.hops.emplace_back ("L","S",-1,1,t_contactL);
        para.hops.emplace_back ("R","S",1,-1,t_contactR);
        auto ampo = get_ampo (glob_basis, sites, para);
        H = toMPO (ampo);
        cout << "MPO dim = " << maxLinkDim(H) << endl;

        // Initialze MPS
        if (Delta != 0.)
        {
            if (basis == "SC")
            {
                psi = get_ground_state_BdG_scatter (glob_basis, sites, mu_biasL, mu_biasS, mu_biasR, para, maxCharge);
            }
            else if (basis == "real_space")
            {
                auto DMRG_sweeps = iut::Read_sweeps (infile, "DMRG_sweeps");
                psi = get_ground_state_SC (glob_basis, sites, mu_biasL, mu_biasS, mu_biasR, para, DMRG_sweeps, args_basis);
            }
            else
            {
                cout << "Unknown basis: " << basis << endl;
                throw;
            }
        }
        else
            psi = get_non_inter_ground_state (glob_basis, sites, mu_biasL, mu_biasS, mu_biasR);
        psi.position(1);

        // Check initial energy
        cout << "Initial energy = " << inner (psi,H,psi) << endl;
    }
    else
    {
        readAll (read_dir+"/"+read_file, psi, H, glob_basis, para, args_basis, step);
        sites = MixedBasis (siteInds(psi), args_basis);
    }
    // ======================= Time evolution ========================
    // Observer
    auto obs = TDVPObserver (sites, psi, {"charge_site",glob_basis.to_glob ("C",1)});
    // Current MPO
    int lenL = visit (basis::size(), glob_basis.parts().at("L"));
    int lenS = visit (basis::size(), glob_basis.parts().at("S"));
    vector<int> spec_links = {lenL-1, lenL+lenS+1};
    int N = glob_basis.N();
    vector<MPO> JMPOs (N);
    for(int i : spec_links)
    {
        auto mpo = get_current_mpo (sites, glob_basis, i);
        JMPOs.at(i) = mpo;
    }
    // Correlation
    /*if (SubCorrN < 0 or SubCorrN > length(psi))
    {
        SubCorrN = length(psi);
    }
    else if (SubCorrN < idevR-idevL+1)
    {
        cout << "sub-correlation size is too small: " << SubCorrN << endl;
        cout << "device size = " << idevL-idevR+1 << endl;
        throw;
    }
    int l = (N - SubCorrN) / 2;
    int ibeg = l+1,
        iend = l+SubCorrN;
    auto sub_corr = SubCorr (glob_basis, ibeg, iend);*/

    // Time evolution
    cout << "Start time evolution" << endl;
    cout << sweeps << endl;
    psi.position(1);
    Real en, err;
    Args args_tdvp_expansion = {"Cutoff",globExpanCutoff, "Method","DensityMatrix",
                                "KrylovOrd",globExpanKrylovDim, "DoNormalize",true, "Quiet",true};
    Args args_tdvp  = {"Quiet",true,"NumCenter",NumCenter,"DoNormalize",true,"Truncate",Truncate,
                       "UseSVD",UseSVD,"SVDmethod",SVDmethod,"WriteDim",WriteDim,"mixNumCenter",mixNumCenter};
    LocalMPO PH (H, args_tdvp);
    while (step <= time_steps)
    {
        cout << "step = " << step << endl;

        // Subspace expansion
        if (maxLinkDim(psi) < sweeps.mindim(1) or (step < globExpanN and (step-1) % globExpanItv == 0))
        {
            timer["glob expan"].start();
            addBasis (psi, H, globExpanHpsiCutoff, globExpanHpsiMaxDim, args_tdvp_expansion);
            PH.reset();
            timer["glob expan"].stop();
        }

        // Time evolution
        timer["tdvp"].start();
        //tdvp (psi, H, 1_i*dt, sweeps, obs, args_tdvp);
        TDVPWorker (psi, PH, 1_i*dt, sweeps, obs, args_tdvp);
        timer["tdvp"].stop();
        auto d1 = maxLinkDim(psi);

        // Measure currents by correlations
        /*timer["current corr"].start();
        sub_corr.measure<MixedBasis> (psi, {"Cutoff",corr_cutoff});
        timer["current corr"].stop();
        for(int j = 1; j < glob_basis.N_phys(); j++)
        {
            timer["current corr"].start();
            auto J = sub_corr.get_current (j);
            timer["current corr"].stop();
            cout << "\t*current " << j << " " << j+1 << " = " << J << endl;
        }*/
        // Measure currents by MPO
        for(int j : spec_links)
        {
            timer["current mps"].start();
            auto J = get_current (JMPOs.at(j), psi);
            timer["current mps"].stop();
            cout << "\t*I " << j << " " << j+1 << " " << J << endl;
        }

        Real NN = den (sites, psi, 1, glob_basis.idevL()-1);
        NN += den (sites, psi, glob_basis.idevL()+1, length(sites));
        NN += den (sites, psi, glob_basis.to_glob ("C",1));
        cout << "tot N = " << NN << endl;

        step++;
        if (write)
        {
            timer["write"].start();
            writeAll (write_dir+"/"+write_file, psi, H, glob_basis, para, args_basis, step);
            timer["write"].stop();
        }
    }
    timer.print();
    return 0;
}
