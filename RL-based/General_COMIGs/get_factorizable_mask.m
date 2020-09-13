function get_factorizable_mask(n, T)

	file = ['Polynomials\polynomialsDeg-3_Terms-', int2str(T), '.mat'];
	
	[~,polys,~] = get_unique_polynomials(n, T);
	factorizable_mask = get_factorizable_polys(polys, n, T);
	save(file, 'factorizable_mask', '-append');
end

function mask = get_factorizable_polys(polys, n, T)

	flag = 0;
	for qubit = 1:n
		flag = flag + check_qubit_factorization(polys, n, T, qubit);
	end
	mask = logical(flag);
end

function mask = check_qubit_factorization(polys, n, T, qubits)
	flag  = ones(size(polys,1),numel(qubits));
	term1 = polys(:, qubits);
	for i = 1:T-1
		flag = flag .* ( term1 == polys(:, qubits + i*n) );
	end
	mask = all(flag,2);
end
