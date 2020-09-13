function Testing_Multiple_Hamiltonians(n, T, c, X)
%	Inputs:
%		n - number of Qubits
%		T - number of Terms
%		c - coefficient constant (coeffs for the quadratic terms in [-c,c])
%		X - optional variable, discard quads that preserve < X states
%

	H = get_unique_polynomials(n,T);
	for i = 1:size(H,1)
		quantum_envelopes_brute_force(H(i,:), c, X);
	end
end
