function out = Find_nullspaces(H)
	[LHS, n] = get_LHS(H);
	[allbits, ~] = get_all_possible_quadratics(n);
	
	allbits_unfolded = reshape( allbits, 2^(2*n), []);
	LHS_allbits_unfolded = reshape( LHS  * allbits, 2^(2*n), []);
	allbits_LHS  = cell2mat( mat2cell( allbits' * LHS, ones(1,size(allbits_unfolded,2))*2^n, 2^n )' );
	allbits_LHS_unfolded = reshape( allbits_LHS, 2^(2*n), []);
	commutator_unfolded = LHS_allbits_unfolded - allbits_LHS_unfolded;
	null_space = null(commutator_unfolded,'r');

	out = size(null_space,2);
end

function [LHS, n] = get_LHS(H, alpha)
	n = max(strlength(H));
	if nargin == 1, alpha = ones(1,size(H,2)); end
	
	LHS = zeros(2^n);
	for t = 1:size(H,2)
		LHS = LHS + alpha(t) * get_term(H{t}, n);
	end
end
