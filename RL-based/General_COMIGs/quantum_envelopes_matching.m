function count = quantum_envelopes_matching(runs, n, LHS, coeffs_all, preserved, allbits_unfolded, terms, null_space, alpha)

[preserved_sorted, I] = sort(preserved(preserved ~= 0),'descend');
coeffs = coeffs_all(I,:);
RHS_cell = mat2cell( reshape( allbits_unfolded * coeffs', 2^n, []), 2^n, ones(1,length(I))*2^n );
[Q, IE] = initialize_sim_diag(LHS);

file_label = 'envelopes/envelopes';
for i = 1:size(alpha,2)
	file_label = [file_label, '_', int2str(alpha(i))];
end
file = fopen([file_label, '.txt'], 'a+');

RHS = cell(1,runs);
count = 0;
for i1 = 1:length(I)
	if preserved_sorted(i1) < 2^n/runs, break, end
	RHS{1} = RHS_cell{i1};
	for i2 = i1+1:length(I)
		if mod(i2, 1000) == 0, fprintf('i1 = %d, i2 = %d\n',i1, i2), end
		idx = [i1, i2];
		if mean(preserved_sorted(idx)) < 2^n/runs, break, end
		RHS{2} = RHS_cell{i2};
		if ~commutes(RHS(1:2)), continue, end
		if runs == 2
			count = get_and_print_envelope(LHS, RHS, runs, Q, IE, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
		else
			for i3 = i2+1:length(I)
				if mod(i3, 1000) == 0, fprintf('i1 = %d, i2 = %d, i3 = %d\n',i1, i2, i3), end
				idx = [i1, i2, i3];
				if mean(preserved_sorted(idx)) < 2^n/runs, break, end
				RHS{3} = RHS_cell{i3};
				if ~commutes(RHS(1:3)), continue, end
				if runs == 3
					count = get_and_print_envelope(LHS, RHS, runs, Q, IE, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
				else
					for i4 = i3+1:length(I)
						idx = [i1, i2, i3, i4];
						if mean(preserved_sorted(idx)) < 2^n/runs, break, end
						RHS{4} = RHS_cell{i4};
						if ~commutes(RHS(1:4)), continue, end
						if runs == 4
							count = get_and_print_envelope(LHS, RHS, runs, Q, IE, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
						else
							for i5 = i4+1:length(I)
								idx = [i1, i2, i3, i4, i5];
								if mean(preserved_sorted(idx)) < 2^n/runs, break, end
								RHS{5} = RHS_cell{i5};
								if ~commutes(RHS(1:5)), continue, end
								if runs == 5
									count = get_and_print_envelope(LHS, RHS, runs, Q, IE, count, coeffs(idx,:), preserved_sorted(idx), terms, null_space, n, file);
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
fprintf('Found %d envelopes!\n', count);
fclose(file);
end
