sigma = cell(4,1);
sigma{1} = [0 1 ; 1 0];
sigma{2} = [0 -1i ; 1i 0];
sigma{3} = [1 0 ; 0 -1];
sigma{4} = eye(2);

x1     = kron(sigma{1},kron(sigma{4},sigma{4}));
x2x3   = kron(sigma{4},kron(sigma{1},sigma{1}));
y2y3   = kron(sigma{4},kron(sigma{2},sigma{2}));
y1z2z3 = kron(sigma{2},kron(sigma{3},sigma{3}));

LHS   = x1*x2x3 + x1*y2y3% + y1z2z3;
QUAD1 =  2*x1 - x2x3 - y2y3 + 2*eye(8)% + y1z2z3;
QUAD2 = -2*x1 + x2x3 + y2y3 + 2*eye(8)% + y1z2z3;

param = initialize_sim_diag(LHS);
param = initialize_sim_diag(param, QUAD1);
[~,D1,D2,D3] = simdiag(param, LHS, QUAD1, QUAD2);

lhs   = real(diag(D1));   assert( all( imag(diag(D1)) < 10e-5 ) );
quad1 = real(diag(D2));   assert( all( imag(diag(D2)) < 10e-5 ) );
quad2 = real(diag(D3));   assert( all( imag(diag(D3)) < 10e-5 ) );

envelope = min(quad1,quad2);
if all( abs(envelope - lhs) < 10e-5 )
	fprintf("Hooray!!!\n");
else
	fprintf(":( :( :(\n");
end
