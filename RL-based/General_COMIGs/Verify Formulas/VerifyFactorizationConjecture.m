n = 10; % number of qubits

sigma = cell(4,1);
sigma{1} = [0 1 ; 1 0];
sigma{2} = [0 -1i ; 1i 0];
sigma{3} = [1 0 ; 0 -1];
sigma{4} = eye(2);

for i=1:3
    for j=1:3
        for k=1:3
			for l=1:3
				for m=1:3
					for p=1:3
						for q=1:3
							for r=1:3
								for s=1:3
									for t=1:3
										LHS = kron(sigma{i},kron(sigma{j},kron(sigma{k},kron(sigma{l},kron(sigma{m},kron(sigma{p},kron(sigma{q},kron(sigma{r},kron(sigma{s},sigma{t})))))))));
										RHS = kron(sigma{i},kron(sigma{j},kron(sigma{k},kron(sigma{l},kron(sigma{m},eye(2^(n-5))))))) - kron(eye(2^5),kron(sigma{p},kron(sigma{q},kron(sigma{r},kron(sigma{s},sigma{t})))));
										RHS1 = eye(2^n) + RHS;
										RHS2 = eye(2^n) - RHS;
										
										[Q, IE] = initialize_sim_diag(LHS);
										[Q, IE] = initialize_sim_diag(Q, IE, RHS1);
										[~,D1,D2,D3] = simdiag(Q, IE, LHS, RHS1, RHS2);
										
										lhs   = real(diag(D1));   assert( all( imag(diag(D1)) < 10e-5 ) );
										quad1 = real(diag(D2));   assert( all( imag(diag(D2)) < 10e-5 ) );
										quad2 = real(diag(D3));   assert( all( imag(diag(D3)) < 10e-5 ) );
										
										envelope  = min(quad1,quad2);
										if all( abs(envelope - lhs) < 10e-5 )
											fprintf("Hooray!!! %d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", i,j,k,l,m,p,q,r,s,t);
										else
											fprintf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d :( :( :(\n", i,j,k,l,m,p,q,r,s,t);
										end
									end
								end
							end
						end
					end
				end
			end
        end
    end
end

