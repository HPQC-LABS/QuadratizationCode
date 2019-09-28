//
//  isquad.cpp
//  

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

// n_max = m_max = 6;
const int P = 4096; // 2^{n+m}
const int Q = 79;   // 1+n+m+\binom{n+m}{2}

// x = 2^j
//
inline int pow2 (int j)
{
  int x = 1;
  for (; j > 0; --j) x *= 2;
  return x;
}

// x_i : where x = (x_0, x_1, ..., x_{n-1})
//       does not ensure that i < n
//
inline int bit (int x, int i)
{
    return (x>>i) & 1;
}

// f(x) = c \prod_{i=1}^n x_i
//
int f (int n, int x, double a[])
{
    int b[n+1];
    for(int i = 1; i <= n; i++){
        b[i] = bit(x, i-1);
    }
    int ans = 0;
    for(int i = 0; i<(1<<n); i++){
        int prod = 1;
        for(int j = 1; j<=n; j++){
            if(((i>>(j-1)) & 1)== 1){
                prod *= b[j];
            }
        }
        ans += a[i] * prod;
    }
    return ans;
}
// ERH : understand why written like this ???
//
string bin2str (int n, int k)
{
  string s;
  for (; n > 0; k = k >> 1, --n) s += '0' + (k % 2);
  return s;
}

//
//
bool isquad (int n, int m, int pn, int pm, double c[], double a[])
{
  int i, j, k, x, y;
  double fx, g, gxy;
  bool tight;
    // bit(x, i-1) = bi
    //c0 + c1 b1 + c2 b2 +... + c4 b4 + c5ba + c6 b1b2 +...
  // ERH : loop over every point
  for (x = 0; x < pn; ++x) {
    fx = f(n,x,a);
    g  = c[0];
    for (i = 1; i <= n; ++i) g += c[i]*(int)bit(x,i-1);

    k = n+m;
    for (i = 1; i < n; ++i)
      for (j = i+1; j <= n; ++j)
        g += c[++k]*(int)bit(x,i-1)*(int)bit(x,j-1);

    // ERH : checking if there is a value of y, so that g(x,y) = f(x)
    tight = false;
    for (y = 0; y < pm; ++y) {
      gxy = g;
      for (j = 1; j <= m; ++j) gxy += c[n+j]*(int)bit(y,j-1);

      k = n+m+n*(n-1)/2;
      for (i = 1; i < m; ++i)
        for (j = i+1; j <= m; ++j)
          gxy += c[++k]*(int)bit(y,i-1)*(int)bit(y,j-1);

      for (i = 1; i <= n; ++i)
        for (j = 1; j <= m; ++j)
          gxy += c[++k]*(int)bit(x,i-1)*(int)bit(y,j-1);

      if (fx == gxy) {
        //cerr << "x = " << bin2str(n,x) << ", y = " << bin2str(m,y) << endl;
        tight = true;
        break;
      }
    }
    if (!tight) return false;
  }
  return true;
}

void print (ostream& ofs, int d, double c[])
{
  for (int i = 0; i < d; ++i)
    ofs << c[i] << " ";
  ofs << endl;
}

double input(){ // new input method on 5 March 2019
    string s;
    cin >> s;
    int l = s.length();
    for(int i = 0; i < l; i++){
        if(s[i] == '/'){
            return stod(s.substr(0, i))/stod(s.substr(i+1, l-i-1));
        }
    }
    return stod(s);
}

int main (int argc, char* argv[])
{
  int d, i, j, m, n;
  char achar;
  double c[Q], a[1<<8];
  ofstream vtx, ray;
  string vtxname, rayname;
  
  if (argc < 2) {
    cout << "isquad <output filename>" << endl;
    return 1;
  }
  vtxname = argv[1]; vtxname += ".gd";
  rayname = argv[1]; rayname += ".bd";
  
  vtx.open (vtxname.c_str());
  ray.open (rayname.c_str());

    cin >> n >> m ;
    for (int k = 0; k < (1<<n); k++){
    cin >> a[k];
}
  d = (1+n+m+((n+m)*(n+m-1)/2));

    vtx << n << " " << m ;
    for (int k = 0; k < (1<<n); k++){
        vtx <<" " << a[k];
    }
    vtx << endl;
    ray << n << " " << m;
    for (int k = 0; k < (1<<n); k++){
        ray <<" " << a[k];
    }
    ray <<endl;
  
  // ERH : for (;;) do this forever 
  for (;;) {
    cin >> achar;
    if (achar == '$') break;
    cin.unget();

    cin >> i; // ignores the first column
    for (i = 0; i < d; ++i) {
        c[i] = input();
    }

    if (isquad (n, m, pow2(n), pow2(m), c, a)) {
      //cerr << "yes" << endl;
      print (vtx, d, c);
    }
    else {
      //cerr << "no" << endl;
      print (ray, d, c);
    }
  }
  vtx << "$" << endl;
  ray << "$" << endl;
  vtx.close();
  ray.close();
  return 0;
}
