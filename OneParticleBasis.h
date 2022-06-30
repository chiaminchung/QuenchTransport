#ifndef __ONEPARTICLEBASIS_H_CMC__
#define __ONEPARTICLEBASIS_H_CMC__
#include "itensor/all.h"
#include "GeneralUtility.h"
using namespace itensor;
using namespace std;

Matrix tight_binding_Hamilt (int L, Real t, Real mu, Real damp_fac=1., bool damp_from_right=true, bool verbose=false)
{
    cout << "L = " << L << endl;
    Matrix H (L,L);
    for(int i = 0; i < L; i++)
    {
        H(i,i) = -mu;
        if (i != L-1)
        {
            int damp_dist = (damp_from_right ? L-2-i : i);
            Real ti = t * pow (damp_fac, damp_dist);
            H(i,i+1) = -ti;
            H(i+1,i) = -ti;
            if (verbose)
                cout << "Hk, t " << i << " = " << ti << endl;
        }
    }
    return H;
}

class OneParticleBasis
{
    public:
        OneParticleBasis () {}
        OneParticleBasis (const string& name, const Matrix& H)
        : _name (name)
        , _H (H)
        {
            diagHermitian (H, _Uik, _ens);
        }
        OneParticleBasis (const string& name, int L, Real t, Real mu, Real damp_fac=1., bool damp_from_right=true, bool verbose=false)
        : _name (name)
        {
            _H = tight_binding_Hamilt (L, t, mu, damp_fac, damp_from_right, verbose);
            diagHermitian (_H, _Uik, _ens);
        }
        OneParticleBasis (const string& name, int L)
        : _name (name)
        {
            _H = Matrix(L,L);
            diagHermitian (_H, _Uik, _ens);
        }

        // Functions that every basis class must have
        const string&                name   ()                const { return _name; }
        vector<tuple<int,auto,bool>> C_op   (int i, bool dag) const;
        Real                         en     (int k)           const { mycheck (k > 0 and k <= _ens.size(), "out of range"); return _ens(k-1); }
        Real                         mu     (int k)           const { mycheck (k > 0 and k <= _ens.size(), "out of range"); return -_H(k-1,k-1); }
        int                          size   ()                const { return _ens.size(); }

        void write (ostream& s) const
        {
            itensor::write(s,_name);
            itensor::write(s,_Uik);
            itensor::write(s,_ens);
        }
        void read (istream& s)
        {
            itensor::read(s,_name);
            itensor::read(s,_Uik);
            itensor::read(s,_ens);
        }

    private:
        string _name;
        Matrix _Uik, _H;
        Vector _ens;
};

// Get the operator information in this basis for the operator Cdag_i, where i is the real-space site index.
// i is 1-index
vector<tuple<int,auto,bool>> OneParticleBasis :: C_op (int i, bool dag) const
{
    mycheck (i > 0 and i <= nrows(_Uik), "out of range");

    // Cdag_i = \sum_k coef_ik Cdag_k
    // Return the basis index (k), coefficient, and whether the operator has a dagger or not
    // k to be returned is 1-index
    auto tmp = _Uik(0,0);
    vector<tuple<int,decltype(tmp),bool>> k_coef_dag;

    for(int k = 0; k < this->size(); k++)       // Here k is zero-index
    {
        k_coef_dag.emplace_back (k+1, iut::conj(_Uik(i-1,k)), dag);
    }
    return k_coef_dag;
}

auto write (ostream& s, const OneParticleBasis& t)
{
    t.write (s);
}
auto read (istream& s, OneParticleBasis& t)
{
    t.read (s);
}
#endif
