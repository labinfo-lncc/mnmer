#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <emscripten.h>

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
            mdat[t.first + " | " + r.first] = r.second;
    

    return mdat;
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

const char* convTab (map<string,string> &mp)
{
    string tab = "<thead><tr><th>mn-mer</th><th>Prob.</th></thead>\n";
    for (auto p: mp)
        tab += "<tr><td>" + p.first + "</td><td>" + p.second + "</td></tr>\n";

    const char* ctab = tab.c_str();

    return ctab;
}

extern "C" {

    const char* EMSCRIPTEN_KEEPALIVE mnmer (const char* seq, int m, int n)
    {
        int k = m + n;
//Get the lexicography of the 4 nucloetide
        vector<string> vp4 = lexnucl(pow(4,k), k);

        map<string,string> mat;

//Header
//        for (auto h: vp4)
//            mat[h] = "0.0";

        map<string,float> mtab = get_mmers (k, m, seq);

        for (auto t: mtab)
            mat[t.first] = to_string(t.second);
        
    
        return convTab(mat);
    }

}

//em++ -std=c++11 mnmercpp.cc -o mnmer.js -O3 -s WASM=1 -s EXPORTED_FUNCTIONS='["_mnmer"]' -s EXPORTED_RUNTIME_METHODS='["cwrap"]'
