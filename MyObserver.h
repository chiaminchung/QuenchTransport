#ifndef __MYOBSERVER_H_CMC__
#define __MYOBSERVER_H_CMC__
#include "itensor/all.h"
#include "Entanglement.h"
using namespace itensor;

template <typename SitesType>
class MyObserver : public DMRGObserver
{
    public:
        MyObserver (const SitesType& sites, const MPS& psi, const Args& args = Args::global())
        : DMRGObserver (psi, args)
        , _sites (sites)
        , _ns (length(psi),0.)
        , _Npar (0.)
        {
            _write = args.getBool ("Write",false);
            _write_minm = args.getInt ("out_minm",0);
            _out_dir = args.getString("out_dir",".");
        }

        void measure (const Args& args = Args::global());

             Real   Npar () const { return _Npar; }
        auto const& ns   () const { return _ns; }

    private:
        bool                 _write;
        string               _out_dir;	// empty string "" if not write
        int                  _write_minm;
        vector<int>          _iDel, _jDel;
        vector<Real>         _Delta;
        SitesType            _sites;

        vector<Real>         _ns;
        Real                 _Npar;
};

Real Onsite_mea (const ITensor& A, const ITensor& op)
{
    ITensor re = A * op;
    re.noPrime ("Site");
    re *= dag(A);
    return re.real();
}

template <typename SitesType>
void MyObserver<SitesType> :: measure (const Args& args)
{
    DMRGObserver::measure (args);

    // Define your measurements below
    // Call psi() to access the MPS
    //
    auto N = length(psi());
    auto b = args.getInt("AtBond",1);
    auto sw = args.getInt("Sweep",0);
    auto ha = args.getInt("HalfSweep",0);
    auto energy = args.getReal("Energy",0);

    // On-site measure
    int oc = orthoCenter(psi());
    if (oc == N || ha == 2) // measure during the second half of sweep
    {
        // Density
        ITensor n_op = _sites.op("N",oc);
        Real ni = Onsite_mea (psi().A(oc), n_op);
        cout << "\tn " << oc << " = " << ni << endl;
        _ns.at(oc-1) = ni;

        // Entanglement entropy
        Real S = EntangEntropy (spectrum());
        cout << "\tentang entropy " << oc << " = " << S << endl;
    }

    if (oc == 1 && ha == 2)
    {
        int m = args.getInt("MaxDim");
        // Write MPS
        if (_write && m >= _write_minm)
        {
            writeToFile (_out_dir+"/psi"+"_m"+to_string(m)+".mps", psi());
        }
    }
}

#endif
