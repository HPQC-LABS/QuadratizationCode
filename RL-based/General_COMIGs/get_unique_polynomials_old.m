function [H, sizes] = get_unique_polynomials(n,T)
%	Inputs:
%		n - number of Qubits
%		T - number of Terms
%	
%	Outputs:
%		H - all unique Hamiltonians with n qubits and T degree-n terms
%		sizes - the sizes of the null spaces of the Hamiltonians in H
%	

file = ['polynomialsDeg-', int2str(n), '_Terms-', int2str(T), '.mat'];
if exist(file,'file'),	load(file, 'H', 'sizes');	return,	end

sz = 3^(n*T);
Hall = 2*ones(1,sz,'uint8');

perm = perms(1:n);
term_perm = cell(1,T);
for i = 1:T, term_perm{i} = perm + (i-1)*n; end
permutations = cell2mat(term_perm(perms(1:T)));

perCheck = 10000;
progress_const = double(floor((sz-1)/perCheck))/100;
progress = 0;
t = cputime;
for checkpoint = 0:floor((sz-1)/perCheck)
	k = int64(checkpoint*perCheck) : min(sz-1, int64((checkpoint+1)*perCheck) - 1 );
	polys = ndec2base(k, 3, n*T);
	if ~mod(checkpoint,10), fprintf('progress = %6.2f%%, time = %6.2f\n', checkpoint/progress_const, cputime - t); end
	for i = k+1
		if Hall(i) == 2
			progress = progress + 1;
			poly = polys(mod(i-1,perCheck)+1,:);
			poly_perms = poly(permutations);
			idx = unique(base2dec(poly_perms,3) + 1);
			Hall(idx)    = 0;
			Hall(idx(1)) = 1;
		end
	end
end
assert(all(Hall ~= 2)); fprintf('progress = %6.2f%%\n',100);
polys = ndec2base(find(Hall)-1, 3, n*T);
polys = remove_duplicate_terms(polys, n, T);
H = get_H(polys, n, T);
save(file, 'polys', 'H');

fprintf('There are %d unique deg-%d polynomials with %d terms\n', size(polys,1), n, T);

%[factors, s1, s2, s3, s4, s5, flag] = find_factors(polys, n, T);
%fprintf('%d of which are factorizable.\n', factors);
%coeffs_not_factorizable = polys(~logical(flag),:);

sizes = get_nullspace_sizes(H, n);
save(file, 'sizes', '-append');
end

function H = get_H(polys, n, T)
%	Splits each poly from a continuous string of Pauli operators
%	into a cell of T terms

	str = char(polys + 40);
	H = cell(size(polys,1),T);
	for poly = 1:size(polys,1)
		for i = 1:T, H{poly,i} = str(poly,(i-1)*n+1:i*n); end
	end
end

function sizes = get_nullspace_sizes(H, n)
%	Prints the mean null space size of the Hamiltonians in H
%
	[allbits, ~] = get_all_possible_quadratics(n);
	just_size = true;
	
	sizes = zeros(1,size(H,1));
	t = cputime;
	for i = 1:size(H,1)
		if ~mod(i,512),	fprintf('progress = %6.2f%%, total time = %4.1f\n',	double(i)/size(H,1)*100, cputime - t);	end
		sizes(i) = find_nullspace(H(i,:), allbits, just_size);
	end
	fprintf('progress = %6.2f%%\n', 100);
	fprintf('avg nullspace = %3.1f, std = %3.1f\n', mean(sizes), std(sizes));
end

function polys = remove_duplicate_terms(polys, n, T)
	idx = zeros(size(polys,1),1);
	for i = 1:T-1
		for j = i+1:T
			idx = idx + all( polys(:, ((i-1)*n+1) : i*n) == polys(:, ((j-1)*n+1) : j*n), 2);
		end
	end
	polys = polys(~logical(idx), :);
end

function [factors, sum1, sum2, sum3, sum4, sum5, flag_] = find_factors(coeffs_all, n, T)
	if T == 1
		factors = 0;sum1=0;sum2=0;sum3=0;sum4=0;sum5=0;
	else
		sum1 = 0;
		flag_ = 0;
		for i = 1:n
			sum1 = sum1 + factorizable(coeffs_all, n, T, i);
			[~, flag] = factorizable(coeffs_all, n, T, i);
			flag_ = flag_ + flag;
		end
		sum2 = 0;
		%flag_ = 0;
		for i = 1:n-1
			for j = i+1:n
				sum2 = sum2 + factorizable(coeffs_all, n, T, [i,j]);
				%[~, flag] = factorizable(coeffs_all, n, T, [i,j]);
				%flag_ = flag_ + flag;
			end
		end
		sum3 = 0;
		if n > 3
			flag_ = 0;
			for i = 1:n-2
				for j = i+1:n-1
					for k = j+1:n
						sum3 = sum3 + factorizable(coeffs_all, n, T, [i,j,k]);
						[~, flag] = factorizable(coeffs_all, n, T, [i,j,k]);
						flag_ = flag_ + flag;
					end
				end
			end
		end
		sum4 = 0;
		if n > 4
			flag_ = 0;
			for i = 1:n-3
				for j = i+1:n-2
					for k = j+1:n-1
						for l = k+1:n
							sum4 = sum4 + factorizable(coeffs_all, n, T, [i,j,k,l]);
							[~, flag] = factorizable(coeffs_all, n, T, [i,j,k,l]);
							flag_ = flag_ + flag;
						end
					end
				end
			end
		end
		sum5 = 0;
		if n > 5
			flag_ = 0;
			for i = 1:n-4
				for j = i+1:n-3
					for k = j+1:n-2
						for l = k+1:n-1
							for m = l+1:n
								sum5 = sum5 + factorizable(coeffs_all, n, T, [i,j,k,l,m]);
								[~, flag] = factorizable(coeffs_all, n, T, [i,j,k,l,m]);
								flag_ = flag_ + flag;
							end
						end
					end
				end
			end
		end
		factors = sum1 - sum2 + sum3 - sum4 + sum5;
	end
end

function [out, flg] = factorizable(coeffs, n, T, qubit)
	flag = ones(size(coeffs,1),numel(qubit));
	for i = 1:T-1
		for j = i+1:T
			flag = flag .* ( coeffs(:, qubit + (i-1)*n ) == coeffs(:, qubit + (j-1)*n) );
		end
	end
	
	flg = ones(size(coeffs,1),1);
	for i = 1:numel(qubit)
		flg = flg .* flag(:,i);
	end
	out = sum(flg);
end

% non-factorizable
function [out, flg] = factorizable2(coeffs, n, T, qubit)
	flag  =  ones(size(coeffs,1),numel(qubit));
	for i = 1:T-1
		for j = i+1:T
			temp = coeffs(:, qubit + (i-1)*n ) ~= coeffs(:, qubit + (j-1)*n);
			flag  = flag .* temp;
		end
	end
	
	flg  = ones(size(coeffs,1),1);
	for i = 1:numel(qubit)
		flg  = flg  .* flag (:,i);
	end
	out = sum(flg);
end

