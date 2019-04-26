//
//  qcone.cpp
//  to run: ./qcone n m a[0] a[1] ... a[2^n - 1] --lrs >temp.ine

#include <cstdio>
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cstring>

using namespace std;

// n_max = m_max = 6;
const int P = 4096; // 2^{n+m}
const int Q = 79;   // 1+n+m+\binom{n+m}{2}

// x_i : where x = (x_0, x_1, ..., x_{n-1})
//       does not ensure that i < n
//
// ERH :
// bitwise operators & and <<, shift 1 i places to the left,
// bitwise & with x. returns the value of bit x_i
//
inline int bit (int x, int i)
{
    return (x>>i) & 1;
}

// x = 2^j
//
inline int pow2 (int j)
{
    int x = 1;
    for (; j > 0; --j) x *= 2;
    return x;
}


// input f here. the bits are b[1..n]
//
// ERH : Calculation of f: x is coded in binary, e.g. when n = 2,
// x = 0 is (0,0), x = 1 is (0, 1), x = 2 is (1, 0), x = 3 is (1, 1)
// comparing x == (pow2(n)-1) is comparing if x == 3, that is (1, 1) and
// therefore, the value of the monomial is 1 too
//
int f (int n, int x, int a[])
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

// ERH : p and q are pn (2^n) and pm (2^m) respectively, in the input
//
void print_fx (int n,  int p, int q, int a[])
{
    int k = 0;
    cout << "param fx := ";
    for (int x = 0; x < p; ++x)
        for (int y = 0; y < q; ++y)
            cout << ++k << " " << f(n, x, a) << " ";
    cout << ";" << endl << endl;
}


// rnd: generates a random integer in [1..x]
//
inline int rnd (int x) { return 1+(int)((double)x*rand()/(RAND_MAX+1.0)); }

// rnd: generates a random integer in [-x..x]
//
inline int srnd (int x)
{
  if (rnd(10)>=5)
    return rnd (x+1)-1;
  else
    return -rnd(x);
}




// matrix form of:
//
// g(x,y) = a + \sum_{i=1}^n b_i x_i
//            + \sum_{j=1}^m c_j y_j
//            + \sum_{1\leq i < j\leq n}  d_{i,j} x_i x_j
//            + \sum_{1\leq i < j\leq m}  e_{i,j} y_i y_j
//            + \sum_{i=1}^n \sum_{j=1}^m f_{i,j} x_i y_j
//
// ERH : in the input, p = pn*pm, q = d
void print_cmatrix (int n, int m, int p, int q, int lrs, int a[])
{
  int i, j, k, l, x, y, p2n, p2m, b[P];
  bool M[P][Q];
  
  p2n = pow2(n);
  p2m = pow2(m);
  
  k = 0;
  // ERH : binary coded x and y -> e.g. n = 2, 
  //       x = 0 -> (0,0)
  //       x = 1 -> (0,1)
  //       x = 2 -> (1,0)
  //       x = 3 -> (1,1)
  // For each of the values that points x, y can take, we substitute 
  // in the equation, some coefficients disappear because x_i or y_j 
  // is 0
    
    // columns are 1 b1 b2 ... bn | ba_1 ... ba_m| b1b2 b1b3 ... b{n-1}bn | ba_1ba_2 ... ba_{m-1}ba_m |b1ba_1 ... bn ba_m
  for (x = 0; x < p2n; ++x)
    for (y = 0; y < p2m; ++y) {
      l = 0;
      M[k][l++] = 1;
      for (i = 0; i < n; ++i) M[k][l++] = bit(x,i);
      for (j = 0; j < m; ++j) M[k][l++] = bit(y,j);
      
      for (i = 0; i < n-1; ++i)
        for (j = i+1; j < n; ++j)
          M[k][l++] = bit(x,i) && bit(x,j);
      
      for (i = 0; i < m-1; ++i)
        for (j = i+1; j < m; ++j)
          M[k][l++] = bit(y,i) && bit(y,j);
      
      for (i = 0; i < n; ++i)
        for (j = 0; j < m; ++j)
          M[k][l++] = bit(x,i) && bit(y,j);
      
      // finally we compute value of f in the point to get the ineq array
        b[k] = f(n,x, a);
      ++k;
    }
  
  // ERH : lrs format -> first independent term (coef of base element "1" - f(x)), then coefs of other elements of the base
  if (lrs) {
    for (i = 0; i < p; ++i) {
      cout << -b[i] << " ";
      for (j = 0; j < q; ++j)
        cout << M[i][j] << " ";
      cout << endl;
    }
  }
  else {
    cout << "param xy: " << endl;
    for (j = 0; j < q; ++j) cout << j << " ";
    cout << ":=" << endl;
    for (i = 0; i < p; ++i) {
      cout << i+1 << " ";
      for (j = 0; j < q; ++j)
        cout << M[i][j] << " ";
      cout << endl;
    }
    cout << ";" << endl << endl;
  }
}

// expression form of:
//
// g(x,y) = a + \sum_{i=1}^n b_i x_i
//            + \sum_{j=1}^m c_j y_j
//            + \sum_{1\leq i < j\leq n}  d_{i,j} x_i x_j
//            + \sum_{1\leq i < j\leq m}  e_{i,j} y_i y_j
//            + \sum_{i=1}^n \sum_{j=1}^m f_{i,j} x_i y_j
//
void print_g (int n, int x, int m, int y)
{
  int i, j, k = 0;

  cout << "c[0]";
  for (i = 0; i < n; ++i) { ++k; if (bit(x,i)) cout << "+c[" << k << "]"; }
  for (j = 0; j < m; ++j) { ++k; if (bit(y,j)) cout << "+c[" << k << "]"; }

  for (i = 0; i < n-1; ++i)
    for (j = i+1; j < n; ++j) {
      ++k; if (bit(x,i) && bit(x,j)) cout << "+c[" << k << "]";
    }

  for (i = 0; i < m-1; ++i)
    for (j = i+1; j < m; ++j) {
      ++k; if (bit(y,i) && bit(y,j)) cout << "+c[" << k << "]";
    }
  
  for (i = 0; i < n; ++i)
    for (j = 0; j < m; ++j) {
      ++k; if (bit(x,i) && bit(y,j)) cout << "+c[" << k << "]";
    }
}

// print_obj
//
void print_obj (int L, int d)
{
  int i;
  int *alpha = new int[d];

  for (i = 0; i < d; ++i) alpha[i] = rnd(L+1)-1;
  for (i = 0; i < d; ++i)
    cout << (alpha[i]==0 ? "+": "") << showpos << alpha[i] << noshowpos << "*c[" << i << "]";

  delete[] alpha;
}

//
//
void print_a (int d, int L)
{
  cout << "param a:= ";
  for (int i = 0; i < d; ++i) cout << i << " " << rnd(L+1)-1 << " ";
  cout << ";" << endl;
}

//
//
void print_aa (int p, int q, int L)
{
  int i, j;
  cout << "param aa : " << endl;
  for (int i = 0; i < q; ++i) cout << i << " ";
  cout << ":=" << endl;
  for (i = 0; i < p; ++i) {
    cout << i+1 << " ";
    for (j = 0; j < q; ++j)
      cout << rnd(L+1)-1 << " ";
    cout << endl;
  }
  cout << ";" << endl << endl;
}

// symmetry breaking
//
void break_symmetry (int n, int m)
{
  int i;

  for (i = 1; i < n; ++i)     cout << "subject to sb" << i << ": c[" << i << "] <= c[" << i+1 << "];" << endl;
  for (i = n+1; i < m+n; ++i) cout << "subject to sc" << i << ": c[" << i << "] <= c[" << i+1 << "];" << endl;
}

void break_symmetry_lrs (int n, int m, int d)
{
  int i, j;
  
  for (i = 1; i < n; ++i) {
    for (j = -1; j < i; ++j) cout << "0 ";
    cout << "-1 1 ";
    for (j = i+2; j < d; ++j) cout << "0 ";
    cout << endl;
  }
    
  for (i = n+1; i < m+n; ++i) {
    for (j = -1; j < i; ++j) cout << "0 ";
    cout << "-1 1 ";
    for (j = i+2; j < d; ++j) cout << "0 ";
    cout << endl;
  }
}

int main (int argc, char* argv[])
{
  // ERH : what is L variable ???
  //       what are preamble commands ??? vble verb, option --nocmd
  int d, k, m, n, x, y, pn, pm, verb = 1, L = 10, hum=0, lrs=0;
  int c;
    int a[1<<8];

  n = 4, m = 1;
  // ERH : initializes random number generator
  srand ((unsigned int) time (0));
  c = rnd (20);
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &m);
    for(int i = 0; i<(1<<n); i++){
        sscanf(argv[i+3], "%d", &a[i]);
    }
  for (--argc; argc > 0; --argc)
    if (!strcmp(argv[argc], "--help")) {
      cout << "qcone : conical polyhedron description of quadratic pBfs" << endl << endl;
      cout << "  --help                 : this message" << endl;
      cout << "  --n=%d           (  4) : number of x variables" << endl;
      cout << "  --m=%d           (  1) : number of y variables" << endl;
      cout << "  --c=%d           (rnd) : coefficient of f(x)" << endl;
      cout << "  --nocmd                : does not print prerambule commands" << endl;
      cout << "  --L=%d           ( 10) : coefficients in [-L..L]" << endl;
      cout << "  --hum            ( nd) : human format" << endl;
      cout << "  --lrs            ( nd) : lrs format" << endl;
      cout << endl;
      cout << "note: if --cpos and --cneg are used simultaneously, only the last one read will take effect" << endl;
      exit (0);
    }
    else if (sscanf(argv[argc], "--n=%d", &n)) assert (n >= 1);
    else if (sscanf(argv[argc], "--m=%d", &m)) assert (m >= 1);
    else if (sscanf(argv[argc], "--c=%d", &c));
    else if (sscanf(argv[argc], "--L=%d", &L));
    else if (!strcmp (argv[argc], "--nocmd")) verb = 0;
    else if (!strcmp (argv[argc], "--hum")) hum = 1;
    else if (!strcmp (argv[argc], "--lrs")) lrs = 1;

  pn = pow2(n);
  pm = pow2(m);
  // ERH : d dimension of the base of polynomials in the expression of g, 
  //       1+n+m+\binom{n+m}{2} 
  d  = (1+n+m+((n+m)*(n+m-1)/2));

  if (lrs) {
    cout << "n" << n << "m" << m << "k" << 3 << endl;
    cout << "H-representation" << endl;
    cout << "begin" << endl;
    cout << pn*pm /*+n+m-2*/ << " " << d+1 << " rational" << endl;
    print_cmatrix (n, m, pn*pm, d, 1, a);
    // break_symmetry_lrs (n, m, d);
    cout << "end" << endl;
    return 0;
  }
  
  if (verb) cout << "option solver cplexamp;" << endl;
  cout << "param DIM := " << d << ";" << endl;
  cout << "param CTR := " << pn*pm << ";" << endl;
  cout << "var c{0..DIM-1};" << endl << endl;

  if (!hum) {
    cout << "param xy{1..CTR,0..DIM-1};" << endl;
    cout << "param fx{1..CTR};" << endl;
    cout << "param a{0..DIM-1};" << endl;
    cout << "param aa{1..CTR,0..DIM-1};" << endl;
    cout << "set XSET;" << endl;
    cout << "data;" << endl;
    print_cmatrix (n, m, pn*pm, d, 0, 0);
    print_fx (n, c, pn, &pm);
    print_a (d, L);
    print_aa (pn*pm, d, L);
    cout << "model;" << endl;

    //    cout << "minimize h{i in 1..CTR}:" << endl;
    //    cout << "  sum{j in 0..DIM-1} xy[i,j]*c[j];" << endl << endl;
    //    cout << "  sum{j in 0..DIM-1} a[j]*c[j];" << endl << endl;
    //    cout << "  sum{j in 0..DIM-1} a[j]*xy[i,j]*c[j];" << endl << endl;
    //    cout << "  sum{j in 0..DIM-1} aa[i,j]*xy[i,j]*c[j];" << endl << endl;
    //    cout << "  sum{j in 0..DIM-1} aa[i,j]*c[j];" << endl << endl;

    cout << "minimize h:" << endl;
    cout << "  sum{i in {1,3,5,9,17,8,12,14,16,20,22,24,26,28,30,32}, j in 0..DIM-1} xy[i,j]*c[j];" << endl << endl;
  }
  else {
    cout << "minimize h:" << endl;
    print_obj (L, d);
    cout << ";" << endl << endl;
  }

  if (!hum) {
    cout << "subject to t{i in 1..CTR}: sum{j in 0..DIM-1} xy[i,j]*c[j] >= fx[i];" << endl;
    
  }
  else {
    k = 0;
    for (x = 0; x < pn; ++x)
      for (y = 0; y < pm; ++y) {
        cout << "subject to t_" << x << "_" << y << "_" << ++k << ": ";
        print_g (n, x, m, y);
        cout << ">=" << f(n, x, a) << ";" << endl;
      }
  }
  // break_symmetry (n, m);
  
  if (verb) {
    //    cout << "solve h" << (hum ? ";" : "[48];") << endl;
    //    cout << "display h" << (hum ? ";" : "[48];") << endl;
    cout << "solve h;" << endl;
    cout << "display h;" << endl;
    cout << "display c;" << endl;
    //    cout << "display c.unbdd;" << endl;
  }
}
