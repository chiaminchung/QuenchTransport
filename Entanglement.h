#ifndef __ENTANGLEMENT_H_CMC__
#define __ENTANGLEMENT_H_CMC__

Real EntangEntropy (const Spectrum& spec)
{
    Real S = 0.;
    for(int i = 1; i <= spec.size(); i++)
    {
        auto p = spec.eig(i);
        S += - p * log(p);
    }
    return S;
}

Real EntangEntropy (const ITensor& C)
{
    auto i1 = C.inds()(1);
    ITensor U(i1), S, V;
    auto spec = svd (C, U, S, V);
    return EntangEntropy (spec);
}
#endif
