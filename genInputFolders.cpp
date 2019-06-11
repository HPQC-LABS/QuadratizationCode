// 5-var with deg 3, 4, 5 terms
// to compile: g++ genInputFolders.cpp -o genInputFolders
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unistd.h>
#include <sys/stat.h>
#include <ctime>

using namespace std;

int pow(int x, int y){
    int ans = 1;
    while(y > 0){
        y--;
        ans *= x;
    }
    return ans;
}

int idtoterm[17];
int termtoid[3000];

vector<int> transform(vector<int> a, int perm[]){
    vector<int> b(a);
    for (int i = 1; i<=10; i++){
        int x[3];
        x[0]= idtoterm[i]/100;
        x[1] = (idtoterm[i]/10) % 10;
        x[2] = idtoterm[i]%10;
        for(int i = 0; i<=2; i++){
            x[i] = perm[x[i]];
        }
        
        sort(x, x+3);
        int term = x[0]*100 + x[1]*10 + x[2];
        b[i] = a[termtoid[term]];
    }
    for (int i = 11; i<=15; i++){
        int x[4];
        x[0]= idtoterm[i]/1000;
        x[1] = (idtoterm[i]/100) % 10;
        x[2] = (idtoterm[i]/10) % 10;
        x[3] = idtoterm[i]%10;
        for(int i = 0; i<=3; i++){
            x[i] = perm[x[i]];
        }
        sort(x, x+4);
        int term = x[0]*1000 + x[1]*100 + x[2]*10 + x[3];
        b[i] = a[termtoid[term]];
    }
    return b;
}

void genInputFolders(int X=500, int N=1, int filehasno=0){
    time_t init, now;
    init = time(NULL);
    idtoterm[1]=123;
    idtoterm[2]=124;
    idtoterm[3]=125;
    idtoterm[4]=134;
    idtoterm[5]=135;
    idtoterm[6]=145;
    idtoterm[7]=234;
    idtoterm[8]=235;
    idtoterm[9]=245;
    idtoterm[10]=345;
    idtoterm[11]=1234;
    idtoterm[12]=1235;
    idtoterm[13]=1245;
    idtoterm[14]=1345;
    idtoterm[15]=2345;
    for (int i = 1; i<=15; i++) {
        termtoid[idtoterm[i]]=i;
    }
    // X the number cores, N the range of coeff, file has no. if filehasno = 1
    set<vector<int> > allcases;
    vector<int> a;

    for (int k = 0; k < pow(2*N+1, 16); k++){
        if (k%10000 == 0) {
            now = time(NULL);
            cerr << k << "  " << allcases.size() << " elapsed time: " <<now - init<<"s"<< endl;
        }
        a.clear();
        a.resize(17);
        //a[0] is not used;
        //a[1..10] are cubic coeff a123, a124, a125, a134, ..., a345;
        //a[11..15] are quartic coeff a1234, a1235, a1245, a1345, a2345
        for (int j = 1; j<=16; j++){
            a[j] = (k / pow(2*N+1, j-1)) % (2*N+1) - N;
        }
        
        //use symmetry
        int perm[] = {0, 1, 2, 3, 4, 5};
        int isnew = 1;
        do {
            vector<int> b = transform(a, perm);
            if(allcases.count(b)!=0){
                isnew = 0;
                break;
            }
        } while ( next_permutation(perm+1,perm+6) );
        int var5 = 1;
        // check whether b1 appear
        if (a[1] == 0 && a[2] == 0 && a[3] == 0 && a[4] == 0 && a[5] == 0 && a[6] == 0 && a[11] == 0 && a[12] == 0 && a[13] == 0 && a[14] == 0){
            var5 = 0;
        }
        // check whether b2 appear
        if (a[1] == 0 && a[2] == 0 && a[3] == 0 && a[7] == 0 && a[8] == 0 && a[9] == 0 && a[11] == 0 && a[12] == 0 && a[13] == 0 && a[15] == 0){
            var5 = 0;
        }
        
        if (a[1] == 0 && a[4] == 0 && a[5] == 0 && a[7] == 0 && a[8] == 0 && a[10] == 0 && a[11] == 0 && a[12] == 0 && a[14] == 0 && a[15] == 0){
            var5 = 0;
        }
        
        if (a[2] == 0 && a[4] == 0 && a[6] == 0 && a[7] == 0 && a[9] == 0 && a[10] == 0 && a[11] == 0 && a[13] == 0 && a[14] == 0 && a[15] == 0){
            var5 = 0;
        }
        
        if (a[3] == 0 && a[5] == 0 && a[6] == 0 && a[8] == 0 && a[9] == 0 && a[10] == 0 && a[12] == 0 && a[13] == 0 && a[14] == 0 && a[15] == 0){
            var5 = 0;
        }
        if (isnew == 1 && var5 == 1){
            allcases.insert(a);
        }
    }
    
    int l = allcases.size();
    
    set<vector<int> >:: iterator it = allcases.begin();
    char temp[50];
    for (int fileno = 1; fileno<= X; fileno++){
        sprintf(temp, "job_%d", fileno);
        mkdir(temp, 0777);
        chdir(temp);
        FILE * f;
        if (filehasno == 1){
            sprintf(temp, "input_coeff_%d.txt", fileno);
            f = fopen(temp, "w");
        }
        else{
            f = fopen("input_coeff.txt", "w");
        }
        for (int k = 1 + (fileno-1) * l/ X; k <= fileno * l/ X; k++){
            a = *it;
            fprintf(f, "0,0,0,0,0,0,0,%d,0,0,0,%d,0,%d,%d,%d,0,0,0,%d,0,%d,%d,%d,0,%d,%d,%d,%d,%d,%d,%d\n", a[1], a[2], a[4], a[7], a[11], a[3], a[5], a[8], a[12], a[6], a[9], a[13], a[10], a[14], a[15], a[16]);
            it++;
        }
        chdir("..");
    }
    return;
}

int main(int argc, char **argv){
    int X, N, filehasno;
    if (argc == 1){
        genInputFolders();
    }
    if (argc == 2){
        sscanf (argv[1], "%d", &X);
        genInputFolders(X);
    }
    if (argc == 3){
        sscanf (argv[1], "%d", &X);
        sscanf (argv[2], "%d", &N);
        genInputFolders(X, N);
    }
    if(argc == 4){
        sscanf (argv[1], "%d", &X);
        sscanf (argv[2], "%d", &N);
        sscanf (argv[3], "%d", &filehasno);
        genInputFolders(X, N, filehasno);
    }
    return 0;
}
