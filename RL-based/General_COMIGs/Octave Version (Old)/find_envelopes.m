function [] = find_envelopes(checkpoint)
	global allbits allbits_unfolded LHS perCheck coeffs_all preserved data_size base coeffs_size n

	k = int64(checkpoint*perCheck) : min( floor(data_size-1), int64((checkpoint+1)*perCheck) - 1 );
	coeffs = ndec2base(k,base,coeffs_size) - '0';
	coeffs( coeffs == 2 ) = -1;
	
	RHS_unfolded = allbits_unfolded * coeffs';
	RHS = sparse( reshape( RHS_unfolded , 2^n, []) );
	
	LHS_RHS = LHS  * RHS;
	RHS_LHS = RHS' * LHS;
	
	LHS_RHS = reshape( permute( reshape(LHS_RHS,2^n,2^n,[]), [1,3,2] ), [], 2^n );
	
	flag = any( abs( LHS_RHS - RHS_LHS ) > 10e-5 , 2);
	for i = 1:n
		flag = flag(1:2:end) + flag(2:2:end);
	end

	if ~all(flag)
		idx = find(flag == 0);
		coeffs_new = coeffs(idx,:);
		coeffs_all{checkpoint + 1} = coeffs_new;

		RHS_new = allbits * kron(coeffs_new',eye(2^n));
		RHS_cell_new = mat2cell(RHS_new,2^n,2^n*ones(1,numel(idx)));
		
		preserved{checkpoint} = zeros(1,numel(idx));
		for i = 1:numel(idx)
			[~,D1,D2] = simdiag(LHS, RHS_cell_new{i});

			lhs = real(diag(D1));  assert( all( imag(diag(D1)) < 10e-5 ) );
			rhs = real(diag(D2));  assert( all( imag(diag(D2)) < 10e-5 ) );

			const = max(lhs - rhs);
			rhs = rhs + const;
			mask_preserve = (abs(rhs - lhs) < 10e-5);
			preserved{checkpoint}(i) = sum(mask_preserve);
		end
	end
end
