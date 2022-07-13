#ifndef __TDVPOBSERVER_H_CMC__
#define __TDVPOBSERVER_H_CMC__
#include <iomanip>
#include <map>
#include "itensor/all.h"
#include "Entanglement.h"
#include "ContainerUtility.h"
using namespace vectool;
using namespace iutility;

template <typename SitesType>
class TDVPObserver : public DMRGObserver
{
    public:
        TDVPObserver (const SitesType& sites, const MPS& psi, const Args& args=Args::global())
        : DMRGObserver (psi, args)
        , _sites (sites)
        , _ns (length(psi),0.)
        , _Npar (0.)
        , _specs (length(psi))
        {
            _write = args.getBool ("Write",false);
            _out_dir = args.getString("out_dir",".");
            _charge_site = args.getInt("charge_site",-1);
        }

        void measure (const Args& args);

             Real   Npar () const { return _Npar; }
        auto const& ns   () const { return _ns; }
        const Spectrum& spec (int i) const { return _specs.at(i); }

    private:
        bool        _write;
        string      _out_dir;	// empty string "" if not write
        SitesType   _sites;
        int         _charge_site=-1;

        // Observables
        vector<Real>        _ns;
        Real                _Npar;
        vector<Spectrum>    _specs;
};

template <typename SitesType>
void TDVPObserver<SitesType> :: measure (const Args& args)
{
    DMRGObserver::measure (args);

    cout << scientific << setprecision(14);
    // Define your measurements below
    // Call psi() to access the MPS
    //
    auto N = length(psi());
    auto b = args.getInt("AtBond");
    auto sw = args.getInt("Sweep");
    auto ha = args.getInt("HalfSweep");
    auto energy = args.getReal("Energy",0);

    if (b != N)
        _specs.at(b) = spectrum();

    int oc = orthoCenter(psi());
    int nc = args.getInt("NumCenter");
    // measure during the second half of sweep
    if (oc == N || ha == 2)
    {
        // Density
        ITensor n_op = noPrime (psi().A(oc) * _sites.op("N",oc), "Site");

        n_op *= dag(psi().A(oc));
        Real ni = real(eltC(n_op));
        cout << "\t*den " << oc << " " << ni << endl;
        _ns.at(oc-1) = ni;

        // Entanglement entropy
        Real S = EntangEntropy (spectrum());
        cout << "\t*entS " << oc << " " << S << endl;

        if (oc == _charge_site)
        {
            auto denmat = psi()(oc) * dag(prime(psi()(oc),"Site"));
            denmat.takeReal();
            auto [Q,D] = diagHermitian (denmat);
            auto ii = denmat.inds()(1);
            int maxOcc = _sites.maxOcc();
            for(int i = 1; i <= dim(ii); i++)
                cout << "\t*nC " << i-maxOcc-1 << " " << elt (denmat,i,i) << endl;
        }
    }

    // At the end of a sweep
    if (oc == 1 && ha == 2 && b == 1)
    {
        for(int i = 1; i < N; i++)
            cout << "\t*m " << i << " " << dim(rightLinkIndex (psi(), i)) << endl;

        if (_write)
        {
            cout << "write MPS" << endl;
            writeToFile (_out_dir+"/psi.mps", psi());
        }
    }
}
#endif
