#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

ifstream in("GCF_002270725.1_ViralProj399942_genomic.fna");
ofstream out("output2.txt");

void back(char c[], string &aux, vector<string> &res, int l) {
    if (l<0)
        res.push_back(aux);
    else {
        for (int i=0; i<4; i++) {
            aux+=c[i];
            back(c, aux, res, l-1);
            aux.pop_back();
        }
    }
}

int main() {

    char c[4]={'A', 'C', 'G', 'T'};
    vector <string> res;
    vector <pair <int, string>> final;
    char s;
    string t;
    int cnt=0;

    while (in >> s) {
        t+=s;
    }

    for (int l=2; l<5; l++) {
        string aux="";
        back(c, aux, res, l);
    }

    for (auto i : res) {

        int ok=0;

        for (int j=0; j+i.size()<=t.size(); j++) {

            string u="";

            for (int k=0; k<i.size(); k++) {
                u+=t[j+k];
            }

            if (u==i)
                ok++;
        }

        if (ok>=6 && ok<=10)
            final.push_back({ok,i});
    }


    for (auto [i, j] : final)
        out<<i<<" "<<j<<endl;

}