n = 4; % number of qubits

sigma = cell(4,1);
sigma{1} = [0 1 ; 1 0];
sigma{2} = [0 -1i ; 1i 0];
sigma{3} = [1 0 ; 0 -1];
sigma{4} = eye(2);

for i=1:4
    for j=1:3
        for k=1:3
			for l=1:3
				LHS  = kron(sigma{i},kron(sigma{j},kron(sigma{k},sigma{l})));
				RHS  = kron(sigma{i},kron(sigma{j},eye(4))) - kron(eye(4),kron(sigma{k},sigma{l}));
				RHS1 = eye(2^n) + RHS;
				RHS2 = eye(2^n) - RHS;
				
				[~,D1,D2,D3] = simdiag(RHS1, RHS2, LHS);
				
				quad1 = real(diag(D1));   assert( all( imag(diag(D1)) < 10e-5 ) );
				quad2 = real(diag(D2));   assert( all( imag(diag(D2)) < 10e-5 ) );
				lhs   = real(diag(D3));   assert( all( imag(diag(D3)) < 10e-5 ) );
				
				envelope  = min(quad1,quad2);
				if all( (envelope - lhs) < 10e-5 )
					fprintf("Hooray!!! %d,%d,%d,%d\n", i,j,k,l);
				else
					fprintf("%d,%d,%d,%d :( :( :(\n", i,j,k,l);
				end
			end
        end
    end
end


