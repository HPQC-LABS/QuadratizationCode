function count = quantum_envelopes_matching(runs, n, LHS, coeffs_all, preserved, allbits_unfolded, terms, null_space, H, c, X)

[preserved, I] = sort(preserved(preserved ~= 0),'descend');
if preserved(1) < 2^n/runs, count = 0; fprintf('Found %d envelopes with %d runs!\n', count, runs);	return;	end
coeffs = coeffs_all(I,:);
I = 1:find(preserved >= X, 1,'last');
if numel(I) < 10000, checkpoint = false; else checkpoint = true; end
RHS_cell = mat2cell( reshape( allbits_unfolded * coeffs(I,:)', 2^n, []), 2^n, ones(1,length(I))*2^n );

file_label = 'envelopes_Partially_Factorizable_Trinomials';
%for i = 1:size(alpha,2), file_label = [file_label, '_', int2str(alpha(i))]; end
file = fopen([file_label, '.txt'], 'a+');
print(H, file);

threshold = 2^n;
sum_mask = zeros(1,10);
hashmap = cell(1,length(I)); imax = 10;
RHS = cell(1,runs); mask = RHS;
count = 0;
param = initialize_sim_diag(LHS);
lhs = param.d;
for i1 = 1:length(I)
	if runs*preserved(i1) < threshold, break, end
	if checkpoint
		if mod(i1, 1) == 0, fprintf('run1 = %d\n',i1), end
	end
	RHS{1} = RHS_cell{i1};
	param1 = initialize_sim_diag(param, RHS{1}); rhs = param1.d;
	rhs = rhs + max(lhs - rhs);
	mask{1}= abs(rhs - lhs) < 10e-5;
	sum_mask(1) = sum(mask{1});
	
	for i2 = i1+1:length(I)
		idx = [i1, i2];
		if sum_mask(1) + (runs-1)*preserved(i2) < threshold, break, end
		if checkpoint
			if mod(i2, 100) == 0, fprintf('run1 = %d, run2 = %d\n',i1, i2), end
		end
		RHS{2} = RHS_cell{i2};
		idx_return = commutes_not(RHS(1:2), hashmap, idx, 1);
		if idx_return, continue, end
		if runs == 2
			count = get_and_print_envelope(LHS, RHS, runs, param1, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
		else
			param2 = initialize_sim_diag(param1, RHS{2}); rhs = param2.d;
			rhs = rhs + max(lhs - rhs);
			mask{2} = logical( mask{1} + ( abs(rhs - lhs) < 10e-5 ) );
			sum_mask(2) = sum(mask{2});
			
			for i3 = i2+1:length(I)
				idx = [i1, i2, i3];
				if sum_mask(2) + (runs-2)*preserved(i3) < threshold, break, end
				if checkpoint
					if ~mod(i3, 10000), fprintf('run1 = %d, run2 = %d, run3 = %d\n',i1, i2, i3), end
				end
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
				if idx_return, continue, end
				if runs == 3
					count = get_and_print_envelope(LHS, RHS, runs, param2, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
				else
					param3 = initialize_sim_diag(param2, RHS{3}); rhs = param3.d;
					rhs = rhs + max(lhs - rhs);
					mask{3} = logical( mask{2} + ( abs(rhs - lhs) < 10e-5 ) );
					sum_mask(3) = sum(mask{3});
					
					for i4 = i3+1:length(I)
						idx = [i1, i2, i3, i4];
						if sum_mask(3) + (runs-3)*preserved(i4) < threshold, break, end
						if checkpoint
							if mod(i4, 10000) == 1, fprintf('run1 = %d, run2 = %d, run3 = %d, run4 = %d\n',i1, i2, i3, i4), end
						end
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
						if idx_return, continue, end
						if runs == 4
							count = get_and_print_envelope(LHS, RHS, runs, param3, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
						else
							param4 = initialize_sim_diag(param3, RHS{4}); rhs = param4.d;
							rhs = rhs + max(lhs - rhs);
							mask{4} = logical( mask{3} + ( abs(rhs - lhs) < 10e-5 ) );
							sum_mask(4) = sum(mask{4});
							
							for i5 = i4+1:length(I)
								idx = [i1, i2, i3, i4, i5];
								if sum_mask(4) + (runs-4)*preserved(i5) < threshold, break, end
								if checkpoint
									if mod(i5, 100000) == 1, fprintf('run1 = %d, run2 = %d, run3 = %d, run4 = %d, run5 = %d\n',i1, i2, i3, i4, i5), end
								end
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
								if idx_return, continue, end
								if runs == 5
									count = get_and_print_envelope(LHS, RHS, runs, param4, count, coeffs(idx,:), preserved(idx), terms, null_space, n, file);
								else
									fprintf('Not implemented for runs > 5\n');
								end
							end
						end
					end
				end
			end
		end
	end
end
fprintf('Found %d envelopes with %d runs!\n', count, runs);
if count
	fprintf(file, 'Found %d envelopes with %d runs! Checked in range [-%d, %d]\n\n', count, runs, c, c);
elseif runs == 3
	fprintf(file, 'Found no envelopes with 2 or 3 runs! Checked in range [-%d, %d]\n\n', c, c);
end
fclose(file);
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

function print(H, file, alpha)
	if nargin == 2, alpha = ones(1,size(H,2)); end
	
	for i = 1:size(H,2)
		if alpha(i) == 1
			fprintf(file, ' + %s', H{i});
		else
			fprintf(file, ' %2d%s', alpha(i), H{i});
		end
	end
	fprintf(file,':\n');
end

