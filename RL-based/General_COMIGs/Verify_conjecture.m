warning('off');
n = 4; % number of qubits

sigma = cell(4,1);
sigma{1} = [0 1 ; 1 0];
sigma{2} = [0 -1i ; 1i 0];
sigma{3} = [1 0 ; 0 -1];
sigma{4} = eye(2);

H = {'XZZY', 'XXZY', 'YXZY', 'YYZY'};
N_of_terms = size(H,2);
terms = cell(1, N_of_terms);
for t = 1:N_of_terms
	terms{t} = get_term(H{t}, n);
end

c = 10;
flag = true;
for i = -c:c
	fprintf('%d\n',i);
    for j = -c:c
        for k = -c:c
			for l = -c:c
				alpha = [i, j, k, l];
				LHS = zeros(2^n);
				for t = 1:N_of_terms
					LHS = LHS + alpha(t) * terms{t};
				end

				coef = sum(abs(alpha));
				RHS  =(coef * eye(2^n) - LHS) * kron(eye(2^(n-2)),kron(sigma{3}, sigma{2}));
				RHS1 = coef * eye(2^n) - RHS;
				RHS2 = coef * eye(2^n) + RHS;

				[Q, IE] = initialize_sim_diag(LHS);
				[~,D1,D2,D3] = simdiag(Q, IE, LHS, RHS1, RHS2);

				lhs   = real(diag(D1));   assert( all( imag(diag(D1)) < 10e-5 ) );
				quad1 = real(diag(D2));   assert( all( imag(diag(D2)) < 10e-5 ) );
				quad2 = real(diag(D3));   assert( all( imag(diag(D3)) < 10e-5 ) );

				envelope  = min(quad1,quad2);
				if ~all( abs(envelope - lhs) < 10e-5 )
					fprintf("%d,%d,%d :( :( :(\n", i,j,k);
					flag = false;
				end
			end
		end
	end
end

if flag, fprintf('HOORAY!!!!\n'), else, fprintf(':( :( :(\n'), end

