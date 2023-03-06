#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <unordered_set>
#include <map>
#include <R.h>
#include <Rdefines.h>

using namespace std;

//Get the number of sequences 
inline unsigned num_seqs (string file)
{
    ifstream f (file);
    string s;
    unsigned n = 0;
    while (getline (f,s))
        if (s[0] == '>')
            ++n;
    f.close();
    return n;
}

inline bool checkseq (string sq, float pni)
{
    if (sq == "") return false;
    
    unsigned sz = sq.size();
    unsigned n = 0;
    unordered_set<char> snuc ( {'A','T','C','G','a','t','c','g'} );

    for (auto c: sq){
        if (snuc.find (c) != snuc.end())
            ++n;
        if ( (float(n)/float(sz)) >= pni )
            return false;
    }

    return true;
    
}

map<string,string> read_fasta_rand (string file, int size, float pni)
{
    unordered_set<unsigned> su;
    mt19937 eng (chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<unsigned> unif (0, size);

    map<string,string> mps;

    while (su.size() != size)
        su.insert (unif(eng));

    unsigned i = 0;
    ifstream f (file);
    string s, id, seq = "";

    while (getline (f,s))
        if (s[0] == '>'){
            
            if ( checkseq (seq,pni) && (su.find(i-1) != su.end()) ){
                mps[id] = seq;
                ++i;              
            }
            
            if (i == size){
                f.close();
                return mps;
            }

            id = s.substr (1,s.find(" ")-1);
            seq = "";
        }       
        else
            seq += s;

    if (su.find(i-1) != su.end())
        mps[id] = seq;

    f.close();

    return mps;
}


extern "C" 
{
    SEXP readrandFASTA (SEXP rfile, SEXP rsize, SEXP rpni)
    {
        string  file = CHAR(STRING_ELT(rfile,0));
        unsigned size = asInteger (rsize);
        float pni = asReal (rpni);

        unsigned num = num_seqs (file);

        if (num < size) return "";

        map<string,string> mpse = read_fasta_rand (file, size, pni);

        string sf = "";

        for (auto p: mpse)
            sf += p.first + "\t" + p.second + "\n";
        
        SEXP cstr = allocVector(STRSXP, sf.size());
        
        PROTECT (cstr);
        cstr = mkString(sf.c_str());
  	    UNPROTECT(1);

        return cstr;
    }
}
