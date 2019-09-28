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

inline int bit(int bitmask, int i){
    return (bitmask>>i) & 1;
}

vector<int> flipbit(vector<int> a, int i){
    vector<int> b(a);
    //a[0] is not used;
    //a[1..10] are cubic coeff a123, a124, a125, a134, ..., a345;
    //a[11..15] are quartic coeff a1234, a1235, a1245, a1345, a2345
    //a[16] coeff a12345
    if (i==1){
        b[16] = -a[16];
        b[15] = a[15] + a[16];
        b[14] = -a[14];
        b[13] = -a[13];
        b[12] = -a[12];
        b[11] = -a[11];
        b[10] = a[10] + a[14];
        b[9] = a[9]+ a[13];
        b[8] = a[8] + a[12];
        b[7] = a[7] + a[11];
        b[6] = -a[6];
        b[5] = -a[5];
        b[4] = -a[4];
        b[3] = -a[3];
        b[2] = -a[2];
        b[1] = -a[1];
    }
    else{
        int perm[] = {0, 1, 2, 3, 4, 5};
        swap(a[1], a[i]);
        a = transform(a, perm);
        a = flipbit(a, 1);
        return transform(a, perm);
    }
    return b;
}

vector<int> flip(vector<int> a, int bitmask){
    //cerr<<"trying flipping "<< bitmask<<endl;
    vector<int> b(a);
    bitmask <<= 1;
    for(int i = 1; i<=5; i++){
        if (bit(bitmask, i) == 1){
            b = flipbit(b, i);
        }
    }
    return b;
}


int main(int argc, char **argv){
    //hardwired data
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
    
    set<vector<int> > allcases;
    vector<int> a;
    a.resize(17);
    //input
    int count = 0;
    while(scanf("0,0,0,0,0,0,0,%d,0,0,0,%d,0,%d,%d,%d,0,0,0,%d,0,%d,%d,%d,0,%d,%d,%d,%d,%d,%d,%d\n", &a[1], &a[2], &a[4], &a[7], &a[11], &a[3], &a[5], &a[8], &a[12], &a[6], &a[9], &a[13], &a[10], &a[14], &a[15], &a[16])==16){
        cerr << allcases.size() <<endl;
        int isnew = 1;
        for (int bitmask = 0; bitmask <= 31; bitmask++){
            if (allcases.count(flip(a, bitmask))!=0){
                isnew = 0;
                cerr<<"already have this if we flip "<< bitmask <<endl;
                break;
            }
        }
        if (isnew == 1){
            allcases.insert(a);
            cerr<<"this is new"<<endl;
            printf("0,0,0,0,0,0,0,%d,0,0,0,%d,0,%d,%d,%d,0,0,0,%d,0,%d,%d,%d,0,%d,%d,%d,%d,%d,%d,%d\n", a[1], a[2], a[4], a[7], a[11], a[3], a[5], a[8], a[12], a[6], a[9], a[13], a[10], a[14], a[15], a[16]);
        }
    }
    return 0;
}
