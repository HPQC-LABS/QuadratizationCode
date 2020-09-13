f_n(D,A,B,C,Y,X) are the formulas for the n=[1,16] coefficients. 

We want to test it for:
-2 -2 -2 -1 0 -1: 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 -6 -6 -7 -8 -5
-2 -2 -1 -1 0 -2: 0 0 0 0 0 0 25 0 0 0 0 0 0 0 0 0 0 -7 -7 -8 -9 -6

We need: 
f_7(-2,-2,-2,-1,0,-1)=24
f_7(-2,-2,-1,-1,0,-2)=25

f_12(-2,-2,-2,-1,0,-1)=-6
f_12(-2,-2,-1,-1,0,-2)=-7

etc.

You pointed out that the formulas are always linear in the coefficients, so each function f_n can be done by matrix-vector multiplication:
[D A B C Y X]*[a1 ; a2 ; a3 ; a4 ; a5 ; a6], where a_m are the 6 coefficients for the function f_1,
[D A B C Y X]*[b1 ; b2 ; b3 ; b4 ; b5 ; b6], where b_m are the 6 coefficients for the function f_2,
etc. up to f_n.

This is then a matrix*matrix multiplication:
LHS_coefficients = NxM where N=# of cases to satisfy, M = # of coefficients in LHS
F = MxL where M = # of coefficients in LHS, and L = # of coefficients in the quad on the RHS.

We would get an NxL matrix, which can be compared to the RHS_coefficients that are desired.

So the formula could be tested for all  ~300 cases with 1 matrix*matrix multiplication and 1 "isequal" application.
For the example above, we want to test:

0 0 0 0 0 0 -3D-3A-3B-3C-3Y-4X 0 0 0 0 0 0 0 0 0 0 , D+B+C+Y+X , D+A+C+Y+X,  D+A+B+Y+X,  D+A+B+C+X,  A+B+C+Y+X

for 
-2 -2 -2 -1 0 -1: 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 -6 -6 -7 -8 -5
-2 -2 -1 -1 0 -2: 0 0 0 0 0 0 25 0 0 0 0 0 0 0 0 0 0 -7 -7 -8 -9 -6

LHS_coefficients = [...
-2 -2 -2 -1 0 -1;
-2 -2 -1 -1 0 -2]

F = [...
0 0 0 0 0 0 -3 0 0 0 0 0 0 0 0 0 0  1 1 1 1 0; % for which coefficients does D appear in the formula, and with what dependence?
 0 0 0 0 0 0 -3 0 0 0 0 0 0 0 0 0 0  0 1 1 1 1; % for which coefficients does A appear in the formula, and with what dependence?
 0 0 0 0 0 0 -3 0 0 0 0 0 0 0 0 0 0  1 0 1 1 1; % for which coefficients does B appear in the formula, and with what dependence?
 0 0 0 0 0 0 -3 0 0 0 0 0 0 0 0 0 0  1 1 0 1 1; % for which coefficients does C appear in the formula, and with what dependence?
 0 0 0 0 0 0 -3 0 0 0 0 0 0 0 0 0 0  1 1 1 0 1; % for which coefficients does Y appear in the formula, and with what dependence?
 0 0 0 0 0 0 -4 0 0 0 0 0 0 0 0 0 0  1 1 1 1 1] % for which coefficients does X appear in the formula, and with what dependence?

Here we get:  LHS_coefficients*F =

     0     0     0     0     0     0    25     0     0     0     0     0     0     0     0     0     0    -6    -6    -7    -8    -6
     0     0     0     0     0     0    26     0     0     0     0     0     0     0     0     0     0    -6    -7    -7    -8    -6


Which can be compared easily to:

0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 -6 -6 -7 -8 -5
0 0 0 0 0 0 25 0 0 0 0 0 0 0 0 0 0 -7 -7 -8 -9 -6

using isequal( LHS_coefficients*F, RHS_coefficients).
