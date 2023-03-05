#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <set>
#include <map>
//#include <R.h>
//#include <Rdefines.h>

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

map<string,string> read_fasta_rand (string file, unsigned num, int mx)
{
    unsigned mx = 100;
    unsigned num = 40;

    if (mx < num) return -1;

    set<unsigned> su;
    mt19937 eng (chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<unsigned> unif (0, mx);

    map<string,string> mps;

    while (su.size() != num)
        su.insert (unif(eng));

    unsigned i = 0;
    ifstream f (file);
    string s, id, seq = "";

    while (getline (f,s))
        if (s[0] == '>'){
            if ( (su.find(i-1) != su.end()) && (seq != "") )
                mps[id] = seq;
            
            id = s.substr (1,s.find(" ")-1);
            seq = "";
            ++i;          
        }       
        else
            seq += s;

    if (su.find(i-1) != su.end())
        mps[id] = seq;


    return 0;
}