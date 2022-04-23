#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

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


void save_file (ofstream &fo, string rname, vector<string> &vp, map<string,float> &mpta)
{
    map<string,float>::iterator it;

    fo <<rname <<flush;
    
    for (auto p: vp)
        if ((it=mpta.find (p)) != mpta.end())
            fo <<',' <<it->second <<flush;
        else
            fo <<",0.0" <<flush;
        fo <<endl;


}


int main (int argc, char* argv[])
{
    if (argc != 7) { cout <<"nmmercpp -f <fasta file> -k <number> -m <number>" <<endl; return 0; }

    int k, m;
    string infile;

    for (int j = 1; argv[j] != NULL; ++j)
        if (string(argv[j]) == "-f")
            infile = argv[++j];
        else if (string(argv[j]) == "-k")
            k = atoi(argv[++j]);
        else if (string(argv[j]) == "-m")
            m = atoi(argv[++j]);


    if (k <= m) { cout <<"Erro in k value. It is equal or greater than m" <<endl; return 0; }

//Get the lexicography of the 4 nucloetide
    vector<string> vp4 = lexnucl(pow(4,k), k);

//Get the sequences
    vector<pair<string,string>> vseqs = get_seqs (infile);

//Result
    ofstream fo ("matrix_nmmer_" + to_string(k-m) + "_" + to_string(m) + ".csv");
//Header
    fo <<"seqid" <<flush;
    for (auto h: vp4)
        fo <<',' <<h <<flush;
    fo <<endl;


#pragma omp parallel for
    for (int i = 0; i < vseqs.size(); ++i){
        map<string,float> mtab = get_mmers (k, m, vseqs[i].second);

#pragma omp critical 
        save_file (fo, vseqs[i].first, vp4, mtab);

    }

    fo.close();

    cout <<"OK" <<endl;

    return 0;
}
