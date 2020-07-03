n = 8; % number of qubits

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
								LHS = kron(sigma{i},kron(sigma{j},kron(sigma{k},kron(sigma{l},kron(sigma{m},kron(sigma{p},kron(sigma{q},sigma{r})))))));
								RHS = kron(sigma{i},kron(sigma{j},eye(2^(n-2)))) - kron(eye(4),kron(sigma{k},kron(sigma{l},kron(sigma{m},kron(sigma{p},kron(sigma{q},sigma{r}))))));
								RHS1 = eye(2^n) + RHS;
								RHS2 = eye(2^n) - RHS;

								[~,D1,D2,D3] = simdiag(RHS1, RHS2, LHS);

								d1 = real(diag(D1));
								d2 = real(diag(D2));
								d3 = real(diag(D3));

								not_preserved = (min(d1,d2) - d3) > 10e-5;
								if ~sum(not_preserved)
									fprintf("Hooray!!! %d,%d,%d,%d,%d,%d,%d,%d\n", i,j,k,l,m,p,q,r);
								else
									fprintf("%d,%d,%d,%d,%d,%d,%d,%d :( :( :(\n", i,j,k,l,m,p,q,r);
								end
							end
						end
					end
				end
			end
        end
    end
end


