#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <R.h>
#include <Rdefines.h>

using namespace std;

inline map<string,float> get_mmers (int k, int m, const string &s)
{
    string q, p1, p2;
    multimap<string,string> mmp;
    set<string> sm;

    for (int i = 0; i < s.size()-k; ++i){
        q = s.substr (i,k);
        p1 = q.substr (0,m);
        p2 = q.substr (m,k);

        mmp.insert ({p1,p2});
        sm.insert (p1);
    }

    pair<multimap<string,string>::iterator,multimap<string,string>::iterator> ret;
    multimap<string,string>::iterator it;

    int u;

    float nc = float(sm.size());

    map<string,map<string,float>> mptab;

    for (auto m: sm){
        ret = mmp.equal_range(m);

        map<string,float> mpf;

        u = distance (ret.first,ret.second);

        for (it=ret.first; it != ret.second; ++it)
            mpf[it->second] += 1.0/(float(u)*nc);
        
        mptab[m] = mpf;
    }

    map<string,float> mdat;

    for (auto t: mptab)
        for (auto r: t.second)
            mdat[t.first+r.first] = r.second;
    

    return mdat;
}

//Read fasta file and save the data in the map structure
vector<pair<string,string>> get_seqs (string file)
{
    ifstream f (file);
    string id, s, seq="";
    vector<pair<string,string>> vps;

    while (getline (f,s))
        if (s[0] == '>'){
            if (seq!="")
                vps.push_back ({id,seq});
            id = s.substr (1,s.size());
            seq="";
        }
        else
            seq += s;

    vps.push_back ({id,seq});

    f.close();

    return vps;
}


//Generates base 4 number  
string conv4 (int num)
{
    string st = to_string(num%4);
    int k = num/4;

    while (k > 3){
        st += to_string (k%4);
        k = k/4;
    }
    st = st + to_string (k%4);
    reverse (st.begin(), st.end());
    return st;
}

//Converte para nucleotideo maiusculo
inline string convNUC (string s)
{
    for (int i = 0; i < s.size(); ++i)
        s[i] = (s[i] == '0') ? 'A' : (s[i] == '1') ? 'T' : (s[i] == '2') ? 'C' : (s[i] == '3') ? 'G' : s[i];
    return s;
}

//Generates a map with nucleotide with base 4 number 
vector<string> lexnucl (int num, int k)
{
    string s;
    vector<string>  vnu;
    int j;

    for (int i = 0; i < num; ++i) {
        s = conv4 (i);

        if (s.size()  < k)
            for (j=s.size(); j < k; ++j)
                s = '0' + s;
    
        vnu.push_back (convNUC(s));    
    }

    return vnu;
}


void save_string (string &res, string rname, vector<string> &vp, map<string,float> &mpta)
{
    map<string,float>::iterator it;

    res += rname;
    
    for (auto p: vp)
        if ((it=mpta.find (p)) != mpta.end())
            res += ',' + to_string(it->second);
        else
            res += ",0.0";
    
    res += "\n";
}



extern "C" {

SEXP cmnmer (SEXP filename, SEXP kk, SEXP mm)
{
    string infile = CHAR(STRING_ELT(filename,0));

    int k = asInteger (kk); 
    int m = asInteger (mm);

//Get the lexicography of the 4 nucloetide
    vector<string> vp4 = lexnucl(pow(4,k), k);

//Get the sequences
    vector<pair<string,string>> vseqs = get_seqs (infile);

//Result
    string result = "seqid";
//Header
    for (auto h: vp4)
        result += ',' + h;
    result += "\n";

    for (int i = 0; i < vseqs.size(); ++i){
        map<string,float> mtab = get_mmers (k, m, vseqs[i].second);

        save_string (result, vseqs[i].first, vp4, mtab);
    }

    SEXP Stab = allocVector(STRSXP, result.size());

    PROTECT (Stab);
    Stab = mkString(result.c_str());
  	UNPROTECT(1);

    return Stab;
}

}
