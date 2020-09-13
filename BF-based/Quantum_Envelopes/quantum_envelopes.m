function quantum_envelopes(n, T, c, X)
%	Inputs:
%		n - number of Qubits
%		T - number of Terms
%		c - coefficient constant (coeffs for the quadratic terms in [-c,c])
%		X - optional variable, discard quads that preserve <= X states
%		

	if nargin == 3, X = 0; end
	file = fopen('progress.txt', 'w');
	
	H = get_unique_polynomials(n, T);
	for i = 1:size(H,1)
		quantum_envelopes_brute_force(H(i,:), c, X, i, file);
	end
end

function quantum_envelopes_brute_force(H, c, X, Hidx, file)
%	Inputs:
%		H - a Hamiltonian to get quadratized
%		c - coefficient constant (coeffs for the quadratic terms in [-c,c])
%		X - discard quads that preserve <= X states
%		

	warning('off');

	[LHS, n, T] = get_LHS(H);
	[allbits, terms] = get_all_possible_quadratics(n);
	param = initialize_sim_diag(LHS);
	lhs = param.d;	rhs = zeros(size(lhs));
	J = param.J;
	K = param.K;

	allbits_unfolded = reshape( allbits, 2^(2*n), []);
	LHS_allbits_unfolded = reshape( LHS  * allbits, 2^(2*n), []);
	allbits_LHS = reshape(permute(reshape( allbits' * LHS, 2^n, [], 2^n), [1,3,2] ), 2^(2*n), []);
	allbits_LHS_unfolded = reshape( allbits_LHS, 2^(2*n), []);
	commutator_unfolded = LHS_allbits_unfolded - allbits_LHS_unfolded;
	null_space = null(commutator_unfolded,'r');

	print(H, Hidx, file);
	%print_nullspace(null_space, terms);
	%null_space = null_space(:,[1,3,4]); % restrict the nullspace

	if isempty(null_space)
		fprintf(file, 'Nothing commutes with the LHS! The Null Space is empty.\n');
		return;
	end
	allbits_unfolded = allbits_unfolded * null_space;
	allbits_size = size(allbits_unfolded,2);

	fprintf(file, '%d Hamiltonians commute with LHS\n',size(null_space,2));
	if size(null_space,2) >= 15, return, end
	%if ~ishermitian(LHS), fprintf("Error! The LHS is not Hermitian!\n"); return, end

	%{
	switch size(null_space,2)
		case 1, c = 15;
		case 2, c = 15;
		case 3, c = 15;
		case 4, c = 10;
		case 5, c = 5;
		case 6, c = 3;
		case 7, c = 2;
		case 8, c = 2;
		otherwise, c = 1;
	end
	%}
	base = size(-c:c,2);
	coeffs_size = allbits_size;

	restartId = 0;
	perCheck = 100000;

	data_size = floor(base^coeffs_size);
	progress_const = 100 * perCheck / data_size;

	coeffs_all = zeros(data_size,coeffs_size);
	preserved  = zeros(1,data_size);

	t = cputime;
	t_init = t;
	for checkpoint = restartId : floor( (data_size-1)/perCheck )
		k = int64(checkpoint*perCheck) : min(data_size-1, int64((checkpoint+1)*perCheck) - 1 );
		coeffs = ndec2base(k,base,coeffs_size) - '0';
		if base > 10, coeffs(coeffs > 9) = coeffs(coeffs > 9) - 7; end	% this is a fix in case base > 10
		%{
		% Putting them in order 0, 1, -1, 2, -2, ...
		for i = 2:2*c
			if mod(i,2)
				coeffs(coeffs == i) = (coeffs(coeffs == i) + 1)/2;
			else
				coeffs(coeffs == i) = -coeffs(coeffs == i)/2;
			end
		end
		%}
		RHS_unfolded = allbits_unfolded * coeffs';
		RHS = reshape( RHS_unfolded , 2^n, []);

		coeffs_all(k + 1,:) = coeffs;

		for i = 1:size(coeffs,1)
			RHS_ = RHS(:, 2^n*(i-1)+1:2^n*i );

			% handle degenerate eigenspaces
			A = param.Q'*RHS_*param.Q;
			for it = 1:numel(J)
				j = J(it); k = K(it);
				rhs(j:k) = real(eig(A(j:k,j:k)));
			end

			rhs = rhs + max(lhs - rhs);
			mask_preserve = (abs(rhs - lhs) < 10e-5);
			preserved(checkpoint * perCheck + i) = sum(mask_preserve);
		end

		fprintf(file,'progress %.2f%%, restart id = %d, step time = %.2f, total time = %.0f, found(%d/%d) = %d, found(%d/%d) = %d, found(%d/%d) = %d, max: %d/%d\n',...
			min((checkpoint+1)*progress_const, 100), checkpoint, cputime - t,cputime - t_init, 2^n/2, 2^n,sum(preserved >= 2^n/2), ceil(2^n/3), 2^n, sum(preserved >= 2^n/3), 2^n/4, 2^n, sum(preserved >= 2^n/4), max(preserved), 2^n);
		t = cputime;
	end

	[preserved, I] = sort(preserved,'descend');
	coeffs = coeffs_all(I,:);
	
	I = 1:find(preserved > X, 1,'last');
	if isempty(I),	fprintf(file, 'No quads that preserve more than %d energies.\n', X);	return;	end
	coeffs = coeffs(I,:);
	preserved = preserved(I);

	RHS_cell = mat2cell( reshape( allbits_unfolded * coeffs', 2^n, []), 2^n, ones(1,length(I))*2^n );
	
	count = 0;
	for runs = 2:4
		if preserved(1) < 2^n/runs, fprintf(file, 'Found 0 envelopes with %d runs!\n', runs);	continue;	end
		[count, quads] = quantum_envelopes_matching(runs, n, LHS, coeffs, preserved, RHS_cell, terms, null_space, file);
		if count
			file_path = ['Quads/quads_', int2str(n),'_', int2str(T),'_', int2str(Hidx), '.mat'];
			save(file_path, 'quads', 'runs');
			break;
		end
	end
	if ~count
		file_path = ['Quads/not_found_', int2str(n), '_', int2str(T), '.txt'];
		print(H, Hidx, fopen(file_path, 'a+'));
	end
end

function term = get_term(str, n)
	sigma = cell(4,1);
	sigma{1} = [0 1 ; 1 0];
	sigma{2} = [0 -1i ; 1i 0];
	sigma{3} = [1 0 ; 0 -1];
	sigma{4} = eye(2);
	
	if str(1) == 'I'
		term = sigma{4};
	else
		term = sigma{str(1) - 87};
	end
	for i = 2:n
		if str(i) == 'I'
			term = kron(term, sigma{4});
		else
			term = kron(term, sigma{str(i) - 87});
		end
	end
end

function [allbits, terms] = get_all_possible_quadratics(n, skip_pauli)
	if nargin < 2
		skip_pauli = 'n';
		file = ['Quadratics/all_possible_quadratics_', int2str(n), '.mat'];
		if exist(file, 'file')
			load(file, 'allbits', 'terms');
			return;
		end
	elseif (skip_pauli ~= 'x') && (skip_pauli ~= 'y') && (skip_pauli ~= 'z')
		fprintf(2,'Error! Wrong input on get_all_possible_quadratics(n, skip_pauli)\n');
	end
	
	sigma = cell(4,1);
	sigma{1} = [0 1 ; 1 0];
	sigma{2} = [0 -1i ; 1i 0];
	sigma{3} = [1 0 ; 0 -1];
	sigma{4} = eye(2);
	
	terms = cell(1,4^n);
	term_idx = 0;
	allbits = zeros(2^n,0);
	
	if n == 3
		for i=1:4
			for j=1:4
				for k=1:4
					allbits = [allbits kron(sigma{i},kron(sigma{j},sigma{k}))];
					term = char([]);
					if i~= 4
						term = [term, char(i+119), '1'];
					end
					if j~= 4
						term = [term, char(j+119), '2'];
					end
					if k~= 4
						term = [term, char(k+119), '3'];
					end
					term_idx = term_idx + 1;
					terms{term_idx} = term;
				end
			end
		end
	elseif n == 4
		for i=1:4
			for j=1:4
				for k=1:4
					for m=1:4
						allbits = [allbits kron(sigma{i},kron(sigma{j},kron(sigma{k},sigma{m})))];
						term = char([]);
						if i~= 4
							term = [term, char(i+119), '1'];
						end
						if j~= 4
							term = [term, char(j+119), '2'];
						end
						if k~= 4
							term = [term, char(k+119), '3'];
						end
						if m~= 4
							term = [term, char(m+119), '4'];
						end
						term_idx = term_idx + 1;
						terms{term_idx} = term;
					end
				end
			end
		end
	elseif n == 5
		for i=1:4
			for j=1:4
				for k=1:4
					for m=1:4
						for p=1:4
							allbits = [allbits kron(sigma{i},kron(sigma{j},kron(sigma{k},kron(sigma{m},sigma{p}))))];
							term = char([]);
							if i~= 4
								term = [term, char(i+119), '1'];
							end
							if j~= 4
								term = [term, char(j+119), '2'];
							end
							if k~= 4
								term = [term, char(k+119), '3'];
							end
							if m~= 4
								term = [term, char(m+119), '4'];
							end
							if p~= 4
								term = [term, char(p+119), '5'];
							end
							term_idx = term_idx + 1;
							terms{term_idx} = term;
						end
					end
				end
			end
		end
	end
	[allbits, terms] = KeepQuadratics(allbits, terms, skip_pauli);
	file = ['all_possible_quadratics_', int2str(n), '.mat'];
	save(file, 'allbits', 'terms');
	
	function [allbits_, terms_] = KeepQuadratics(allbits, terms, skip_pauli)
	width = size(allbits,1);
	
	terms_ = cell(1,2); % Size unknown, it will fill up dynamically
	allbits_  = [];
	count = 0;
	for k = 1:size(terms,2)
		switch size(terms{k},2)
			case 2,	flag = (terms{k}(1) ~= skip_pauli);
			case 4, flag = (terms{k}(1) ~= skip_pauli) && (terms{k}(3) ~= skip_pauli);
			otherwise, flag = false;
		end
		if flag
			count = count + 1;
			terms_{count} = terms{k};
			allbits_ = [allbits_ allbits( :, (k-1)*width+1 : k*width ) ];
		end
	end
end
end

function print(H, Hidx, file)

	for i = 1:size(H,2)
		if i ~= 1
			fprintf(file, ' + %s', H{i});
		else
			fprintf(file, '%s', H{i});
		end
	end
	fprintf(file, ', id = %d\n', Hidx);
end

function param = initialize_sim_diag(varargin)
	if nargin == 1
		A = varargin{1};
		Q = eye(size(A));
		IE = [ones(1,length(Q)-1),0];
		J = find((IE==1).*[1,IE(1:end-1)==0]);
		K = [];
		for j = J, K = [K, j+find(IE(j:end)==0,1,'first')-1 ]; end
		epsilon = 10*max(size(A))*eps(normest(A));
		cutoff  = -round(log10(epsilon))-3;
	else
		Q  = varargin{1}.Q;
		IE = varargin{1}.IE;
		J  = varargin{1}.J;
		K  = varargin{1}.K;
		epsilon = varargin{1}.epsilon;
		cutoff  = varargin{1}.cutoff;
		A = varargin{2};
		A = Q'*A*Q;
	end
	
	for i = 1:numel(J)
		j = J(i); k = K(i);
		[V,d] = schur(A(j:k,j:k)); d = diag(d);
		[~,IS] = sort(round(d,cutoff)); % handle case d = [1, i, (1-eps)*i]
		d = d(IS).';
		Q(:,j:k) = Q(:,j:k)*V(:,IS);
		IE(j:k) = [abs(d(1:end-1)-d(2:end)) < epsilon,0];
	end
	param.Q  = Q;
	param.IE = IE;
	param.J  = find((IE==1).*[1,IE(1:end-1)==0]);
	param.K  = [];
	for j = param.J, param.K = [param.K, j+find(IE(j:end)==0,1,'first')-1 ]; end
	param.epsilon = epsilon;
	param.cutoff  = cutoff;
	if nargin == 1
		param.d = real(diag(Q'*A*Q));
	else
		param.d = real(diag(Q'*varargin{2}*Q));
	end
end

function [count, quads] = quantum_envelopes_matching(runs, n, LHS, coeffs, preserved, RHS_cell, terms, null_space, file)

	threshold = 2^n;
	data_size = numel(preserved);

	sum_mask = zeros(1,10);
	hashmap = cell(1,data_size);
	imax = 10;
	RHS = cell(1,runs);
	mask = RHS;
	count = 0;

	param = initialize_sim_diag(LHS);
	lhs = param.d;

	progress = 'run1 = %d';
	for k = 2:runs
		progress = [progress, ', run', int2str(k), ' = %d'];
	end
	progress = [progress, '\n'];

	counter = 0;
	print_threshold = 100000;
	quads = [];
	for i1 = 1:data_size
		if runs*preserved(i1) < threshold, break, end
		if counter > print_threshold,	print_progress(file, progress, runs, i1); counter = 0;	end
		RHS{1} = RHS_cell{i1};
		param1 = initialize_sim_diag(param, RHS{1}); rhs = param1.d;
		rhs = rhs + max(lhs - rhs);
		mask{1}= abs(rhs - lhs) < 10e-5;
		sum_mask(1) = sum(mask{1});

		for i2 = i1+1:data_size
			idx = [i1, i2];
			if sum_mask(1) + (runs-1)*preserved(i2) < threshold, break, end
			if counter > print_threshold,	print_progress(file, progress, runs, idx); counter = 0;	end
			RHS{2} = RHS_cell{i2};
			idx_return = commutes_not(RHS(1:2), hashmap, idx, 1);
			counter = counter + 1;
			if idx_return, continue, end
			if runs == 2
				[count, quad] = get_and_print_envelope(LHS, RHS, runs, param1, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
				quads = [quads; quad];
			else
				param2 = initialize_sim_diag(param1, RHS{2}); rhs = param2.d;
				rhs = rhs + max(lhs - rhs);
				mask{2} = logical( mask{1} + ( abs(rhs - lhs) < 10e-5 ) );
				sum_mask(2) = sum(mask{2});

				for i3 = i2+1:data_size
					idx = [i1, i2, i3];
					if sum_mask(2) + (runs-2)*preserved(i3) < threshold, break, end
					if counter > print_threshold,	print_progress(file, progress, runs, idx); counter = 0;	end
					RHS{3} = RHS_cell{i3};
					idx_return = commutes_not(RHS(1:3), hashmap, idx, 2);
					if idx_return > 0
						if isempty(hashmap{idx(end)})
							hashmap{idx(end)} = zeros(1,imax);
							hashmap{idx(end)}(1) = 1; % counter
							hashmap{idx(end)}(2) = idx(idx_return);
						elseif hashmap{idx(end)}(end) ~= 0
							hashmap{idx(end)}(2*end) = 0;
							hashmap{idx(end)}(1) = hashmap{idx(end)}(1) + 1;
							hashmap{idx(end)}(end/2+1) = idx(idx_return);
						else
							hashmap{idx(end)}(1) = hashmap{idx(end)}(1) + 1;
							hashmap{idx(end)}(hashmap{idx(end)}(1)+1) = idx(idx_return);
						end
					end
					counter = counter + 1;
					if idx_return, continue, end
					if runs == 3
						[count, quad] = get_and_print_envelope(LHS, RHS, runs, param2, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
						quads = [quads; quad];
					else
						param3 = initialize_sim_diag(param2, RHS{3}); rhs = param3.d;
						rhs = rhs + max(lhs - rhs);
						mask{3} = logical( mask{2} + ( abs(rhs - lhs) < 10e-5 ) );
						sum_mask(3) = sum(mask{3});

						for i4 = i3+1:data_size
							idx = [i1, i2, i3, i4];
							if sum_mask(3) + (runs-3)*preserved(i4) < threshold, break, end
							if counter > print_threshold,	print_progress(file, progress, runs, idx); counter = 0;	end
							RHS{4} = RHS_cell{i4};
							idx_return = commutes_not(RHS(1:4), hashmap, idx, 3);
							if idx_return > 0
								if isempty(hashmap{idx(end)})
									hashmap{idx(end)} = zeros(1,imax);
									hashmap{idx(end)}(1) = 1; % counter
									hashmap{idx(end)}(2) = idx(idx_return);
								elseif hashmap{idx(end)}(end) ~= 0
									hashmap{idx(end)}(2*end) = 0;
									hashmap{idx(end)}(1) = hashmap{idx(end)}(1) + 1;
									hashmap{idx(end)}(end/2+1) = idx(idx_return);
								else
									hashmap{idx(end)}(1) = hashmap{idx(end)}(1) + 1;
									hashmap{idx(end)}(hashmap{idx(end)}(1)+1) = idx(idx_return);
								end
							end
							counter = counter + 1;
							if idx_return, continue, end
							if runs == 4
								[count, quad] = get_and_print_envelope(LHS, RHS, runs, param3, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
								quads = [quads; quad];
							else
								param4 = initialize_sim_diag(param3, RHS{4}); rhs = param4.d;
								rhs = rhs + max(lhs - rhs);
								mask{4} = logical( mask{3} + ( abs(rhs - lhs) < 10e-5 ) );
								sum_mask(4) = sum(mask{4});

								for i5 = i4+1:data_size
									idx = [i1, i2, i3, i4, i5];
									if sum_mask(4) + (runs-4)*preserved(i5) < threshold, break, end
									if counter > print_threshold,	print_progress(file, progress, runs, idx); counter = 0;	end
									RHS{5} = RHS_cell{i5};
									idx_return = commutes_not(RHS(1:5), hashmap, idx, 4);
									if idx_return > 0
										if isempty(hashmap{idx(end)})
											hashmap{idx(end)} = zeros(1,imax);
											hashmap{idx(end)}(1) = 1; % counter
											hashmap{idx(end)}(2) = idx(idx_return);
										elseif hashmap{idx(end)}(end) ~= 0
											hashmap{idx(end)}(2*end) = 0;
											hashmap{idx(end)}(1) = hashmap{idx(end)}(1) + 1;
											hashmap{idx(end)}(end/2+1) = idx(idx_return);
										else
											hashmap{idx(end)}(1) = hashmap{idx(end)}(1) + 1;
											hashmap{idx(end)}(hashmap{idx(end)}(1)+1) = idx(idx_return);
										end
									end
									counter = counter + 1;
									if idx_return, continue, end
									if runs == 5
										[count, quad] = get_and_print_envelope(LHS, RHS, runs, param4, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
										quads = [quads; quad];
									else
										fprintf(file, 'Not implemented for runs > 5\n');
									end
								end
							end
						end
					end
				end
			end
		end
	end
	fprintf(file, 'Found %d envelopes with %d runs!\n', count, runs);
	
	function print_progress(file, progress, runs, idx)

		for i = numel(idx)+1:runs
			idx = [idx, idx(end) + 1];
		end
		fprintf(file, progress, idx);
	end
end

function flag = commutes_not(RHS, hashmap, idx, n)
% AB is Hermitian for Hermitian A,B iff A and B commutes
if (n~=1 && idx(1)~=1 && idx(2)~=2)
	if ~isempty(hashmap{idx(end)})
		for i = n:-1:1
			if any(hashmap{idx(end)}(2:end) == idx(i))
				flag = -1;
				return;
			end
		end
	end
end

for i = n:-1:1
	if ~ishermitian( RHS{i} * RHS{end} )
		flag = i;
		return;
	end
end
flag = 0;
end

function [count, quad] = get_and_print_envelope(LHS, RHS, runs, param, count, coeffs, preserved, terms, null_space, n, file)

	[flag, const] = get_envelope;
	if flag
		print_envelope(coeffs, preserved, int8(const));
		quad = [coeffs, const];
		count = count + 1;
	else
		quad = [];
	end
	
	function [flag, const] = get_envelope
		D = cell(1,runs+1);
		d = cell(1,runs+1);
		const = zeros(1,runs);

		[~, D{:}] = simdiag(param, LHS, RHS{:});
		for i = 1:runs+1
			d{i} = real(D{i});
			assert( all( imag(D{i}) < 10e-5 ) );
		end

		for i = 1:runs
			const(i) = max(d{1}-d{i+1});
			d{i+1} = d{i+1} + const(i);
		end

		for i = 3:runs+1
			d{2} = min(d{2},d{i});
		end

		flag = all( (d{2} - d{1}) < 10e-5 );
	end
	function [] = print_envelope(quads, pres, const)
		fprintf(file,"LHS = min(");
		for k = 1:size(quads,1)
			print_quad(quads(k,:), const(k));
			if k ~= size(quads,1)
				fprintf(file,',');
			else
				fprintf(file,' )');
			end
		end
		for k = 1:size(quads,1)
			fprintf(file, ', %d/%d', pres(k), 2^n);
		end
		fprintf(file,'\n');
		
		function [] = print_quad(quad, const)
			for i = 1:numel(quad)
				for j = 1:size(terms,2)
					coef = quad(i) * null_space(j,i);
					if coef == 1
						fprintf(file, ' + %s',terms{j});
					elseif coef == -1
						fprintf(file, ' - %s',terms{j});
					elseif coef > 0
						fprintf(file, ' + %d%s',coef, terms{j});
					elseif coef < 0
						fprintf(file, ' - %d%s',-coef, terms{j});
					end
				end
			end
			if const > 0
					fprintf(file, ' + %d',const);
				elseif const < 0
					fprintf(file, ' - %d',-const);
			end
		end
	end
end

function S = ndec2base(D,B,N)
%NDEC2BASE Convert decimal integer to base B string.
%   NDEC2BASE(D,B) returns the representation of D as a string in
%   base B.  D must be a non-negative integer array smaller than 2^52
%   and B must be an integer between 2 and 36.
%
%   FDEC2BASE(D,B,N) produces a representation with at least N digits.
%
%   NDEC2BASE is designed to be a direct replacement for dec2base. It is
%   optimized for speed and conservation of memory.
% 
%   Examples
%       fdec2base(23,3) returns '212'
%       fdec2base(23,3,5) returns '00212'
%

narginchk(2,3);

D = D(:);
if ~(isnumeric(D) || ischar(D))
    error('First argument must be an array of integers, 0 <= D <= 2^52.');
end
if ~isscalar(B) || ~(isnumeric(B) || ischar(B)) || B ~= floor(B) || B < 2 || B > 36
    error('Second argument must be an integer, 2 <= B <= 36.');
end
if nargin == 3
    if ~isscalar(N) || ~(isnumeric(N) || ischar(N)) || N ~= floor(N) || N < 0
        error('Third argument must be a positive integer.');
    end
end

D = double(D);
B = double(B);
l = length(D);

maxd = max(D);
no = max(1,round(log2(maxd+1)/log2(B)));
while B^no <= maxd, no=no+1; end

if nargin == 3 && N>no, nmax = N;else, nmax = no; end
S = zeros(l,no,'uint8');

% for large arrays a lot faster than matrix multiplication
while no > 1
    no  = no - 1;
    base = B^no;
    col  = nmax-no;
    S(:,col) = uint8(D/base-0.5);
    D  = D - double(S(:,col))*base;
end

S(:,nmax) = uint8(D);
if B>10
    symbols = uint8('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ');
    % conserve memory + faster than s = symbols(s+1);
    S=S+1;
    maxn = size(S,1);
    for a = 1:10000:maxn
        ind = a:min(a+9999,maxn);
        S(ind,:) = symbols(S(ind,:));
    end
    S = char(S);
else
    S = S+uint8('0');
    S = char(S);
end
end

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
		if nargout  > 1,	polys = char(char(join(H,'',2)) - 40);	end
		if nargout == 3,	load(file, 'sizes');	end
		return;
	end

	if T == 1
		polys = ndec2base(0:3^n-1, 3, n);
		polys = remove_duplicate_entries(polys, 1, 3);

		H = get_H(polys + 40, n, T);
		save(file, 'H');
		sizes = get_nullspace_sizes(H, n);
		save(file, 'sizes', '-append');
		factorizable_mask = get_factorizable_polys(polys, n, T);
		save(file, 'factorizable_mask', '-append');
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
		if checkpoint, fprintf('progress = %6.2f%%, time = %6.2f\n', 100*double(checkpoint)/floor(size(polys,1)/perCheck), cputime - t); end

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
	polys = polys(Hall,:);
	fprintf('There are %d unique deg-%d polynomials with %d terms\n', size(polys,1), n, T);

	H = get_H(polys + 40, n, T);
	save(file, 'H');

	sizes = get_nullspace_sizes(H, n);
	save(file, 'sizes', '-append');
	factorizable_mask = get_factorizable_polys(polys, n, T);
	save(file, 'factorizable_mask', '-append');

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
	%
	
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
				for k=1:numelA
					lia(k) = any(a(k)==b(:));
				end
			else
				for k=1:numelA
					found = a(k)==b(:);
					if any(found)
						lia(k) = true;
						locb(k) = find(found,1);
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
					for k=1:numelA
						lia(k) = any(a(k)==b(:));   % ANY returns logical.
					end
				else
					for k=1:numelA
						found = a(k)==b(:); % FIND returns indices for LOCB.
						if any(found)
							lia(k) = true;
							found = find(found);
							locb(k) = found(1);
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

function [LHS, n, T] = get_LHS(H, terms)

	if nargin == 1
		n = max(strlength(H));
		T = size(H,2);
		
		LHS = get_term(H{1}, n);
		for t = 2:T,	LHS = LHS + get_term(H{t}, n);	end
	else
		LHS = terms(H{1});
		for t = 2:size(H,2),	LHS = LHS + terms(H{t});	end
	end
end

function [] = print_nullspace(nullspace, terms)
	
	for i = 1:size(nullspace, 2)
		fprintf('\n%d)', i);
		for j = 1:size(nullspace, 1)
			coef = nullspace(j,i);
			if coef == 1
				fprintf(' + %s',terms{j});
			elseif coef == -1
				fprintf(' - %s',terms{j});
			elseif coef > 0
				fprintf(' + %d%s',coef, terms{j});
			elseif coef < 0
				fprintf(' - %d%s',-coef, terms{j});
			end
		end
	end
	fprintf('\n');	
end

function [Q,varargout] = simdiag(param, varargin)
%
% [Q,D1,...,Dm] = simdiag(A1,...,Am,options);
%	Simultaneously diagonalize all input matrices
%
%	Input:		A1,...,A2	pairwise commuting complex normal matrices
%				options		options field:
%								tol: tolerance
%
%	Ouput:		Q			unitary matrix containing the simultaneous
%							eigenvectors of A1,...,Am
%				D1,...,Dm	eigenvalues of A1,...,Am, respectively
%
%	Copyright (c) 2008-2009, Christian Mendl
%	All rights reserved.

tol = eps^1.5;
n = size(varargin{1},1);

% preprocessing step
Q = dodo(param.Q, param.J, param.K, param.cutoff, varargin{end});
for j=1:length(varargin)
	varargin{j} = Q'*varargin{j}*Q;
end

% Reference:
%	Angelika Bunse-Gerstnert, Ralph Byers, and Volker Mehrmann,
%	Numerical Methods for Simultaneous Diagonalization,
%	SIAM J. Matrix Anal. Appl. Vol. 14, No. 4, pp. 927-949, October 1993
calc_off2 = @(A)norm(A-diag(diag(A)),'fro')^2;
max_iter = 100;
num_iter = 0;
off2 = 0; nscale = 0;
for m=1:length(varargin)
	off2 = off2 + calc_off2(varargin{m});
	nscale = nscale + norm(varargin{m},'fro');
end
while off2 > tol*nscale
	for j=1:n
		for k=j+1:n
			v = zeros(length(varargin),3);
			for m=1:length(varargin)
				v(m,:) = [varargin{m}(j,j)-varargin{m}(k,k),varargin{m}(j,k),varargin{m}(k,j)];
			end
			[c,s] = approx_min(v);
			Q = timesR(Q,j,k,c,s);
			for m=1:length(varargin)
				varargin{m} = rotate(varargin{m},j,k,c,s);
			end
		end
	end
	off2 = 0; nscale = 0;
	for m=1:length(varargin)
		off2 = off2 + calc_off2(varargin{m});
		nscale = nscale + norm(varargin{m},'fro');
	end
	% number of iterations
	num_iter = num_iter + 1;
	if num_iter > max_iter
		fprintf('Exiting: Maximum number of iterations exceeded. Current relative error: %g.\n',off2/nscale);
		break;
	end
end

% eigenvalues
for j=1:min(length(varargin),nargout-1)
	varargout{j} = diag(varargin{j});
end
end

%% A = A*R, R = R(j,k,c,s)
function A = timesR(A,j,k,c,s)

% A = A*R
A(:,[j,k]) = [c*A(:,j)+s*A(:,k),-conj(s)*A(:,j)+conj(c)*A(:,k)];
end

%% A = R'*A*R, R = R(j,k,c,s)
function A = rotate(A,j,k,c,s)

A(:,[j,k]) = [c*A(:,j)+s*A(:,k),-conj(s)*A(:,j)+conj(c)*A(:,k)];	% A = A*R
A([j,k],:) = [conj(c)*A(j,:)+conj(s)*A(k,:);-s*A(j,:)+c*A(k,:)];	% A = R'*A
end

%%
% Approximate minimizer of
% |s c conj(v(:,1)) - c^2 conj(v(:,2)) + s^2 conj(v(:,3))|^2 + |s c v(:,1) + s^2 v(:,2) - c^2 v(:,3)|^2
% for c = cos(theta), s = exp(i phi) sin(theta)
function [c,s] = approx_min(v)

target = @(c,s,v) norm(s*c*conj(v(:,1))-c^2*conj(v(:,2))+s^2*conj(v(:,3)),2)^2+norm(s*c*v(:,1)+s^2*v(:,2)-c^2*v(:,3),2)^2;

[c,s] = calc_min(v(1,1),v(1,2),v(1,3));
m = target(c,s,v);
for j=2:size(v,1)
	[c1,s1] = calc_min(v(j,1),v(j,2),v(j,3));
	x = target(c1,s1,v);
	if x < m
		m = x;
		c = c1; s = s1;
	end
end
end

%%
% Exact minimizer of
% |s c conj(a0) - c^2 conj(a21) + s^2 conj(a12)|^2 + |s c a0 + s^2 a21 - c^2 a12|^2;
% Refer to
%	H. H. Goldstine and L. P. Horwitz, A Procedure for the
%	Diagonalization of Normal Matrices, J. ACM (1959)
function [c,s] = calc_min(a0,a21,a12)

u = real(a0);
v = imag(a0);

tmp = (a21+conj(a12))/2;
r = abs(tmp); beta = angle(tmp);

tmp = (a21-conj(a12))/2;
s = abs(tmp); gamma = angle(tmp);

nu = beta-gamma;
sin_nu = sin(nu);
cos_nu = cos(nu);

L = u*v-4*r*s*sin_nu;
M = u^2-v^2+4*(r^2-s^2);

A = L*(r^2-s^2)*sin_nu+M*r*s;
B = L*(r^2+s^2)*cos_nu;
C = L*(r^2-s^2)+M*r*s*sin_nu;

tmp = r*s*cos_nu*sqrt(M^2+4*L^2);
phi = (atan2(-A*C+B*tmp,B*C+A*tmp)-beta-gamma)/2;

r_cos_ba = r*cos(beta+phi);
s_sin_ga = s*sin(gamma+phi);
kappa = u^2+v^2-4*(r_cos_ba^2+s_sin_ga^2);
lambda = 4*(u*r_cos_ba+v*s_sin_ga);
theta = -atan2(-lambda,kappa)/4;

c = cos(theta);
s = exp(1i*phi)*sin(theta);
end

%% "DODO" (diagonalize one then diagonalize the other) preprocessing step
% note: this is unstable in general!
function Q = dodo(Q, J, K, cutoff, A)
%	handle degenerate eigenspaces

	A = Q'*A*Q;
	for i = 1:numel(J)
		j = J(i); k = K(i);
		[V,d] = schur(A(j:k,j:k));
		[~,IS] = sort(round(diag(d),cutoff)); % handle case d = [1, i, (1-eps)*i]
		Q(:,j:k) = Q(:,j:k)*V(:,IS);
	end
end
