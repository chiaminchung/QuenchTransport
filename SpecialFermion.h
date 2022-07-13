#ifndef __SPECIALFERMION_H_CMC__
#define __SPECIALFERMION_H_CMC__
#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"
using namespace itensor;

class SpecialFermionSite
{
    Index s;
    public:

    SpecialFermionSite(Index I) : s(I) { }

    SpecialFermionSite(Args const& args = Args::global())
    {
        auto systype = args.getString("SystemType");
        auto ts = TagSet("Site,Fermion");
        auto n = 1;
        if(args.defined("SiteNumber"))
        {
          n = args.getInt("SiteNumber");
          ts.addTags("n="+str(n));
        }
        auto in_scatter = args.getBool("in_scatter");
        if (systype == "SC_scatter")
        {
            if (in_scatter)
            {
                s = Index(QN({"Nf",0,-1},{"Ps",0,-2}),1,
                          QN({"Nf",0,-1},{"Ps",1,-2}),1,
                          Out,ts);
            }
            else
            {
                s = Index(QN({"Nf",0,-1},{"Ps",0,-2}),1,
                          QN({"Nf",1,-1},{"Ps",0,-2}),1,
                          Out,ts);
            }
        }
        else if (systype == "Normal")
        {
            if (in_scatter)
            {
                s = Index(QN({"Nf",0,-1},{"Ns",0,-1}),1,
                          QN({"Nf",1,-1},{"Ns",1,-1}),1,
                          Out,ts);
            }
            else
            {
                s = Index(QN({"Nf",0,-1},{"Ns",0,-1}),1,
                          QN({"Nf",1,-1},{"Ns",0,-1}),1,
                          Out,ts);
            }
        }
        else if (systype == "SC_Josephson_scatter")
        {
            if (in_scatter)
            {
                s = Index(QN({"Pf",0,-2},{"Ps",0,-2}),1,
                          QN({"Pf",0,-2},{"Ps",1,-2}),1,
                          Out,ts);
            }
            else
            {
                s = Index(QN({"Pf",0,-2},{"Ps",0,-2}),1,
                          QN({"Pf",1,-2},{"Ps",0,-2}),1,
                          Out,ts);
            }
        }
        else
        {
            cout << "Unknown system type: " << systype << endl;
            throw;
        }
    }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "Emp" || state == "0") 
            {
            return s(1);
            }
        else 
        if(state == "Occ" || state == "1") 
            {
            return s(2);
            }
        else
            {
            throw ITError("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args) const
        {
        auto sP = prime(s);

        auto Emp  = s(1);
        auto EmpP = sP(1);
        auto Occ  = s(2);
        auto OccP = sP(2);
         
        auto Op = ITensor(dag(s),sP);

        if(opname == "N" || opname == "n")
            {
            Op.set(Occ,OccP,1);
            }
        else
        if(opname == "C")
            {
            Op.set(Occ,EmpP,1);
            }
        else
        if(opname == "Cdag")
            {
            Op.set(Emp,OccP,1);
            }
        else
        if(opname == "A")
            {
            Op.set(Occ,EmpP,1);
            }
        else
        if(opname == "Adag")
            {
            Op.set(Emp,OccP,1);
            }
        else
        if(opname == "F" || opname == "FermiPhase")
            {
            Op.set(Emp,EmpP,1);
            Op.set(Occ,OccP,-1);
            }
        else
        if(opname == "projEmp")
            {
            Op.set(Emp,EmpP,1);
            }
        else
        if(opname == "projOcc")
            {
            Op.set(Occ,OccP,1); 
            }
        else
        if(opname == "I")
            {
            Op.set(Occ,OccP,1);
            Op.set(Emp,EmpP,1);
            }
        else
            {
            throw ITError("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
};

using SpecialFermion = BasicSiteSet<SpecialFermionSite>;
#endif
