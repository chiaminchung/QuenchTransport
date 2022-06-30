#ifndef __SPECIALBOSON_H_CMC__
#define __SPECIALBOSON_H_CMC__
#include "itensor/mps/siteset.h"

class SpecialBosonSite
{
    Index s;

    vector<int> _ns;

    public:

    int n (int i) const { return _ns.at(i-1); }

    SpecialBosonSite(Index I, Args const& args = Args::global()) : s(I)
    {
        auto maxOcc = args.getInt("MaxOcc");
        for(int n = -maxOcc; n <= maxOcc; n++)
            _ns.push_back (n);
    }

    SpecialBosonSite(Args const& args = Args::global())
    {
        auto systype = args.getString("SystemType");

        auto tags = TagSet("Site,Boson");
        if (args.defined("SiteNumber"))
        {
            auto n = args.getInt("SiteNumber");
            tags.addTags("n="+str(n));
        }

        auto maxOcc = args.getInt("MaxOcc");
        for(int n = -maxOcc; n <= maxOcc; n++)
            _ns.push_back (n);

        auto qints = Index::qnstorage(_ns.size());
        for(int i = 0; i < _ns.size(); i++)
        {
            int n = _ns.at(i);
            if(systype == "SC_scatter")
            {
                int p = n % 2;
                qints[i] = QNInt(QN({"Nf",n,-1},{"Ps",p,-2}),1);
            }
            else if (systype == "Normal")
            {
                qints[i] = QNInt(QN({"Nf",0,-1},{"Ns",-n,-1}),1);
            }
            else if (systype == "SC_Josephson_scatter")
            {
                int p = n % 2;
                qints[i] = QNInt(QN({"Pf",p,-2},{"Ps",p,-2}),1);
            }
            else
            {
                cout << "Unknown system type: " << systype << endl;
                throw;
            }
        }
        s = Index(std::move(qints),Out,tags);
    }

    Index
    index() const { return s; }

    IndexVal
    state(std::string state)
    {
	    if(state == "Emp")
            state = "0";
        for(int i = 0; i < _ns.size(); i++)
        {
            if(state == str(_ns.at(i))) return s(i+1);
        }
        throw ITError("State " + state + " not recognized");
        return IndexVal{};
    }

	ITensor op (std::string const& opname, Args const& args) const
    {
        auto sP = prime(s);

        auto Op = ITensor(dag(s),sP);

        if(opname == "N" || opname == "n")
        {
            for(int i = 0; i < _ns.size(); i++)
            {
                int j = i+1;
                int n = _ns.at(i);
                Op.set(s=j,sP=j,n);
            }
        }
        else if(opname == "NSqr" || opname == "nSqr")
        {
            for(int i = 0; i < _ns.size(); i++)
            {
                int j = i+1;
                int n = _ns.at(i);
                Op.set(s=j,sP=j,n*n);
            }
        }
        else if(opname == "A" or opname == "C")
        {
            for(int i = 1; i < _ns.size(); i++)
            {
                Op.set(s=1+i,sP=i,1);
            }
        }
        else if(opname == "Adag" or opname == "Cdag")
        {
            for(int i = 1; i < _ns.size(); i++)
            {
                Op.set(s=i,sP=1+i,1);
            }
        }
        else if(opname == "A2")
        {
            for(int i = 3; i <= _ns.size(); i++)
            {
                Op.set(s=i,sP=i-2,1);
            }
        }
        else if(opname == "A2dag")
        {
            for(int i = 3; i <= _ns.size(); i++)
            {
                Op.set(s=i-2,sP=i,1);
            }
        }
        else if(opname == "I")
        {
            for(int i = 1; i <= _ns.size(); i++)
            {
                Op.set(s=i,sP=i,1);
            }
        }
        else
        {
            throw ITError("Operator \"" + opname + "\" name not recognized");
        }

        return Op;
    }
};
#endif
