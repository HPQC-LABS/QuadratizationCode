function count = quantum_envelopes_matching(runs, n, LHS, coeffs_all, preserved, allbits_unfolded, terms, null_space, H, c)

[preserved_sorted, I] = sort(preserved(preserved ~= 0),'descend');
coeffs = coeffs_all(I,:); I = 1:100;
RHS_cell = mat2cell( reshape( allbits_unfolded * coeffs(I,:)', 2^n, []), 2^n, ones(1,length(I))*2^n );

file_label = 'envelopes_monomial';
%for i = 1:size(alpha,2)
%	file_label = [file_label, '_', int2str(alpha(i))];
%end
file = fopen([file_label, '.txt'], 'a+');
%fprintf(file, '%s + %s + %s:\n', H{1}, H{2}, H{3});

threshold = 2^n/runs;
RHS = cell(1,runs);
count = 0;
[Q, IE] = initialize_sim_diag(LHS);
for i1 = 1:length(I)
	if preserved_sorted(i1) < threshold, break, end
	RHS{1} = RHS_cell{i1};
	[Q1,IE1] = initialize_sim_diag(Q, IE, RHS{1});
	for i2 = i1+1:length(I)
		if mod(i2, 1000) == 0, fprintf('i1 = %d, i2 = %d\n',i1, i2), end
		idx = [i1, i2];
		if mean(preserved_sorted(idx)) < threshold, break, end
		RHS{2} = RHS_cell{i2};
		if ~commutes(RHS(1:2)), continue, end
		if runs == 2
			count = get_and_print_envelope(LHS, RHS, runs, Q1, IE1, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
		else
			[Q2,IE2] = initialize_sim_diag(Q1, IE1, RHS{2});
			for i3 = i2+1:length(I)
				if mod(i3, 1000) == 0, fprintf('i1 = %d, i2 = %d, i3 = %d\n',i1, i2, i3), end
				idx = [i1, i2, i3];
				if mean(preserved_sorted(idx)) < threshold, break, end
				RHS{3} = RHS_cell{i3};
				if ~commutes(RHS(1:3)), continue, end
				if runs == 3
					count = get_and_print_envelope(LHS, RHS, runs, Q2, IE2, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
				else
					[Q3,IE3] = initialize_sim_diag(Q2, IE2, RHS{3});
					for i4 = i3+1:length(I)
						idx = [i1, i2, i3, i4];
						if mean(preserved_sorted(idx)) < threshold, break, end
						RHS{4} = RHS_cell{i4};
						if ~commutes(RHS(1:4)), continue, end
						if runs == 4
							count = get_and_print_envelope(LHS, RHS, runs, Q3, IE3, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
						else
							[Q4,IE4] = initialize_sim_diag(Q3, IE3, RHS{4});
							for i5 = i4+1:length(I)
								idx = [i1, i2, i3, i4, i5];
								if mean(preserved_sorted(idx)) < threshold, break, end
								RHS{5} = RHS_cell{i5};
								if ~commutes(RHS(1:5)), continue, end
								if runs == 5
									count = get_and_print_envelope(LHS, RHS, runs, Q4, IE4, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
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
if (count == 0) && (runs == 3)
	fprintf(file, 'Found no envelopes with 2 or 3 runs! Checked in range [-%d, %d]\n\n', c, c);
elseif count
	fprintf(file, 'Found %d envelopes with %d runs! Checked in range [-%d, %d]\n\n', count, runs, c, c);
end
fclose(file);
end

function flag = commutes(RHS)
	flag = all (cellfun(@(x) all( all( abs( RHS{end}*x - x*RHS{end} ) < 10e-5 ) ), RHS(1:end-1)) );
end
