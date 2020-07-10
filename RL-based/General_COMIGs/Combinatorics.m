n = 3; % # of Qubits
T = 3; % # of Terms

sz = 3^(n*T);
Hall = -1*ones(1,sz);
k = 0:sz - 1;
coeffs = ndec2base(k,3,n*T);
perm = perms(1:n);

switch T
	case 1, permutations = perm;
	case 2, permutations = [ [perm, n + perm]; [perm + n, perm] ];
	case 3, permutations = [ [perm, n + perm, 2*n + perm]; [perm, 2*n + perm, n + perm]; [n + perm, perm, 2*n + perm]; [n + perm, 2*n + perm, perm]; [2*n + perm, n + perm, perm]; [2*n + perm, perm, n + perm] ];
end

for i = 1:sz
	if ~mod(i,10000), fprintf('progress = %6.2f%%\n', double(i)/sz*100 ); end
	if Hall(i) == -1
		coef = coeffs(i,:);
		tern = coef(permutations);
		idx = unique(base2dec(tern,3) + 1);
		Hall(idx(1)) = 1;
		for k = 2:numel(idx)
			Hall(idx(k)) = 0;
		end
	end
end
assert(all(Hall ~= -1)); fprintf('progress = %6.2f%%\n',100);
coeffs_all = coeffs(logical(Hall),:);
idx = get_idx(coeffs_all, n, T);
coeffsfinal = coeffs_all(idx,:);

[factors, s1, s2, s3, s4, s5, flag] = find_factors(coeffsfinal, n, T);
fprintf('There are %d unique deg-%d polynomials with %d terms', size(coeffsfinal,1), n, T);
%fprintf(', %d of which are factorizable.\n', factors);

%str = char(coeffsfinal+40);
str = char(coeffsfinal(logical(flag),:)+40);
sizes = [];
for i = 1:size(str,1)
	if mod(i,10) == 1, fprintf('i = %d\n', i); end
	switch T
		case 1, sizes = [sizes, Find_nullspaces( { str(i,1:n) } ) ];
		case 2, sizes = [sizes, Find_nullspaces( { str(i,1:n), str(i,n+1:2*n) } ) ];
		case 3, sizes = [sizes, Find_nullspaces( { str(i,1:n), str(i,n+1:2*n), str(i,2*n+1:3*n) } ) ];
		case 4, sizes = [sizes, Find_nullspaces( { str(i,1:n), str(i,n+1:2*n), str(i,2*n+1:3*n), str(i,3*n+1:4*n) } ) ];
	end
end
fprintf('avg nullspace = %3.1f, std = %3.1f\n', mean(sizes), std(sizes));


function idx_bool = get_idx(coef, n, T)
	idx = zeros(size(coef,1),1);
	for i = 1:T-1
		for j = i+1:T
			idx = idx + all( coef(:, ((i-1)*n+1) : i*n) == coef(:, ((j-1)*n+1) : j*n), 2);
		end
	end
	idx_bool = ~logical(idx);
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

