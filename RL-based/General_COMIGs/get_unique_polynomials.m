function [H, polys, sizes] = get_unique_polynomials(n,T)
%	Inputs:
%		n - number of Qubits
%		T - number of Terms
%		
%	Outputs:
%		H - all unique Hamiltonians with n qubits and T degree-n terms
%		polys - a form of the polynomials used for the recursion
%		sizes - the sizes of the null spaces of the Hamiltonians in H
%		

file = ['Polynomials/polynomialsDeg-', int2str(n), '_Terms-', int2str(T), '.mat'];
if exist(file,'file')
	load(file, 'H');
	if nargout  > 1,	polys = char(char(join(H,'')) - 40);	end
	if nargout == 3,	load(file, 'sizes');	end
	return;
end

[~,polys,~] = get_unique_polynomials(n, T-1);

new_term = ndec2base(0:3^n-1, 3, n);
polys_new = zeros(size(new_term,1) * size(polys,1), n*T);
for i = 1:size(new_term,1)
	polys_new( (i-1)*size(polys,1) + 1 : i*size(polys,1), :) = [polys, repmat(new_term(i,:),size(polys,1),1)];
end
polys = remove_duplicate_terms(polys_new, n, T);	clear polys_new
polys = remove_duplicate_entries(polys, n, T);

perm = perms(1:n);
perm(factorial(n),:) = []; % get rid of the identity perm: 1,2,...,n
term_perm = cell(1,T);
for i = 1:T, term_perm{i} = perm + (i-1)*n; end
permutations = cell2mat(term_perm);

perCheck = 10000;
rows	= size(permutations,1);
columns = size(permutations,2);

poly_perms_all = zeros(perCheck * rows, columns);
Hall = true(1,size(polys,1));
t = cputime;
for checkpoint = 0:floor(size(polys,1)/perCheck)
	for i = checkpoint*perCheck+1:min( (checkpoint+1)*perCheck, size(polys,1))
		poly = polys(i,:);
		poly_perms = poly(permutations);
		poly_perms_all( mod(i-1,perCheck)*rows+1:(mod(i-1,perCheck) + 1)*rows,:) = poly_perms;
	end
	poly_perms_sorted = sort_terms(poly_perms_all, n, T);
	if ~mod(i,1), fprintf('progress = %6.2f%%, time = %6.2f\n', 100*double(checkpoint)/floor(size(polys,1)/perCheck), cputime - t); end
	
	[lia,loc] = ismember_fast(poly_perms_sorted, polys(Hall,:));
	temp = find(Hall);
	for i = checkpoint*perCheck+1:min( (checkpoint+1)*perCheck, size(polys,1))
		if Hall(i)
			lia_ = lia( mod(i-1,perCheck)*rows+1:(mod(i-1,perCheck) + 1)*rows );
			if any(lia_)
				loc_ = loc( mod(i-1,perCheck)*rows+1:(mod(i-1,perCheck) + 1)*rows);
				loc_ = loc_(lia_);
				Hall(temp(loc_)) = 0;
				Hall(i) = 1;
			end
		end
	end
end
fprintf('progress = %6.2f%%, time = %6.2f\n', 100, cputime - t);
fprintf('There are %d unique deg-%d polynomials with %d terms\n', size(polys,1), n, T);

polys = polys(Hall,:);
H = get_H(polys + 40, n, T);
save(file, 'H');

sizes = get_nullspace_sizes(H, n);
save(file, 'sizes', '-append');
factorizable_mask = get_factorizable_polys(polys, n, T);
save(file, 'factorizable_mask', '-append');
end

function [lia, locb] = ismember_fast(a,b)
%	Faster version of 'ismember'
%
	if size(a,1) == 1
		uA = repmat(a, size(b,1), 1);
		d = uA(1:end,:) == b(1:end,:);
		d = all(d,2);
		lia = any(d);
		if lia
			locb = find(d, 1, 'first');
		end
		return;
	end
	
	[uA,~,icA] = unique(a,'rows','sorted');
	[sortuAuB,IndSortuAuB] = sortrows([uA;b]);
	
	% Find matching entries
	d = sortuAuB(1:end-1,:)==sortuAuB(2:end,:);		% d indicates matching entries
	d = all(d,2);                                   % Finds the index of matching entries
	ndx1 = IndSortuAuB(d);                          % NDX1 are locations of repeats in C
	
	[lia, locb] = ismemberBuiltinTypes(icA, ndx1);
	d = find(d);
	newd = d(locb(lia));
	where = IndSortuAuB(newd+1)-size(uA,1);
	locb(lia) = where;
end

function [lia,locb] = ismemberBuiltinTypes(a,b)
% General handling.
% Use FIND method for very small sizes of the input vector to avoid SORT.
if nargout > 1
    locb = zeros(size(a));
end
% Handle empty arrays and scalars.  
numelA = numel(a);
numelB = numel(b);
if numelA == 0 || numelB <= 1
    if numelA > 0 && numelB == 1
        lia = (a == b);
        if nargout > 1
            % Use DOUBLE to convert logical "1" index to double "1" index.
            locb = double(lia);
        end
    else
        lia = false(size(a));
    end
    return
end

scalarcut = 5;
if numelA <= scalarcut
    lia = false(size(a));
    if nargout <= 1
        for i=1:numelA
            lia(i) = any(a(i)==b(:));
        end
    else
        for i=1:numelA
            found = a(i)==b(:);
            if any(found)
                lia(i) = true;
                locb(i) = find(found,1);
            end
        end
    end
else
    % Use method which sorts list, then performs binary search.
    % Convert to full to work in C helper.
    if issparse(a)
        a = full(a);
    end
    if issparse(b)
        b = full(b);
    end
    
    if (isreal(b))
        % Find out whether list is presorted before sort
        sortedlist = issorted(b(:));
        if nargout > 1
            if ~sortedlist
                [b,idx] = sort(b(:));
            end
        elseif ~sortedlist
            b = sort(b(:));
        end
    else
        sortedlist = 0;
        [~,idx] = sort(real(b(:)));
        b = b(idx);
    end
    
    % Use builtin helper function ISMEMBERHELPER:
    % [LIA,LOCB] = ISMEMBERHELPER(A,B) Returns logical array LIA indicating
    % which elements of A occur in B and a double array LOCB with the
    % locations of the elements of A occurring in B. If multiple instances
    % occur, the first occurrence is returned. B must be already sorted.
    
    if ~isobject(a) && ~isobject(b) && (isnumeric(a) || ischar(a) || islogical(a))
        if (isnan(b(end)))
            % If NaNs detected, remove NaNs from B.
            b = b(~isnan(b(:)));
        end
        if nargout <= 1
            lia = builtin('_ismemberhelper',a,b);
        else
            [lia, locb] = builtin('_ismemberhelper',a,b);
        end
    else % a,b, are some other class like gpuArray, sym object.
        lia = false(size(a));
        if nargout <= 1
            for i=1:numelA
                lia(i) = any(a(i)==b(:));   % ANY returns logical.
            end
        else
            for i=1:numelA
                found = a(i)==b(:); % FIND returns indices for LOCB.
                if any(found)
                    lia(i) = true;
                    found = find(found);
                    locb(i) = found(1);
                end
            end
        end
    end
    if nargout > 1 && ~sortedlist
        % Re-reference locb to original list if it was unsorted
        locb(lia) = idx(locb(lia));
    end
end
end

function polys = remove_duplicate_entries(polys, n, T)
	
	polys = sort_terms(polys, n, T);
	polys = unique(polys, 'rows');
end

function polys = sort_terms(polys, n, T)

	polys_splitted = get_H(polys, n, T);
	polys_sorted = sort(polys_splitted,2);
	polys = char(join(polys_sorted,''));
end

function H = get_H(polys, n, T)
%	Splits each poly from a continuous string of Pauli operators
%	into a cell of T terms
	
	c = char(polys);
	H = strings(size(c,1), T);
	for i = 1:T
		idx	= (i-1)*n+1:i*n;
		H(:,i) = string(c(:,idx));
	end
end

function sizes = get_nullspace_sizes(H, n)
%	Prints the mean null space size of the Hamiltonians in H
%
	[allbits, ~] = get_all_possible_quadratics(n);
	terms = get_terms(n);
	
	sizes = zeros(1,size(H,1),'uint8');
	t = cputime;
	for i = 1:size(H,1)
		if ~mod(i,2048),	fprintf('progress = %6.2f%%, time = %4.1f\n',	double(i)/size(H,1)*100, cputime - t);	end
		sizes(i) = find_nullspace(H(i,:), allbits, terms, n);
	end
	fprintf('progress = %6.2f%%, time = %4.1f\n', 100, cputime - t);
	fprintf('avg nullspace = %3.1f, std = %3.1f\n', mean(double(sizes)), std(double(sizes)));
end

function terms = get_terms(n)

	sigma = cell(4,1);
	sigma{1} = [0 1 ; 1 0];
	sigma{2} = [0 -1i ; 1i 0];
	sigma{3} = [1 0 ; 0 -1];
	sigma{4} = eye(2);
	
	k = ndec2base(0:3^n-1, 3, n);
	str_all = k - 47;
	values = cell(1,3^n);
	for s = 1:3^n
		str = str_all(s,:);
		term = sigma{str(1)};
		for i = 2:n,	term = kron(term, sigma{str(i)});	end
		values{s} = term;
	end
	
	keys = string(char(k + 40));
	terms = containers.Map(keys, values);
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

