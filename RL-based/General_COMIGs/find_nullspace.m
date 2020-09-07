function out = find_nullspace(H, allbits, terms, n)
	
	if nargin == 1
		[LHS, n] = get_LHS(H);
		allbits	 = get_all_possible_quadratics(n);
	else
		LHS = get_LHS(H, terms);
	end
	
	LHS_allbits_unfolded = reshape( LHS  * allbits, 2^(2*n), []);
	allbits_LHS = reshape(permute(reshape( allbits' * LHS, 2^n, [], 2^n), [1,3,2] ), 2^(2*n), []);
	allbits_LHS_unfolded = reshape( allbits_LHS, 2^(2*n), []);
	commutator_unfolded = LHS_allbits_unfolded - allbits_LHS_unfolded;
	if nargin == 4
		out = size(commutator_unfolded,2) - rank(commutator_unfolded);
	else
		out = null(commutator_unfolded,'r');
	end
end

function [LHS, n] = get_LHS(H, terms)

	if nargin == 1
		n = max(strlength(H));
		
		LHS = get_term(H{1}, n);
		for t = 2:size(H,2),	LHS = LHS + get_term(H{t}, n);	end
	else
		LHS = terms(H{1});
		for t = 2:size(H,2),	LHS = LHS + terms(H{t});	end
	end
end
