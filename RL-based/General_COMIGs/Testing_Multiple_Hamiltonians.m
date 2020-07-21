n = 3;
T = 2;

coeffs = get_unique_polynomials(n,T);
for i = 1:size(coeffs,1)
	H = get_H(coeffs(i,:), n, T);
	quantum_envelopes_brute_force(H);
end

function H = get_H(coef, n, T)

	str = char(coef+40);
	H = cell(1,T);
	for i = 1:T, H{i} = str((i-1)*n+1:i*n); end
end
