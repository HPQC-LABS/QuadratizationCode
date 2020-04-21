runs = 2;

[preserved_sorted, I] = sort(preserved(preserved ~= -1),'descend');
coeffs = coeffs_all(I,:);

RHS = allbits * kron(coeffs', speye(2^n));
RHS_cell = mat2cell(RHS, 2^n, ones(1,length(I))*2^n);

RHS = cell(1,runs  );
D   = cell(1,runs+1);
d   = cell(1,runs+1);

const = zeros(1,runs);
count = 0;

for k = 1:length(I)-1
	if preserved_sorted(k) < 2^n/2
		break
	end
	RHS{1} = RHS_cell{k};
	for j = k+1:length(I)
		if preserved_sorted(k) + preserved_sorted(j) < 2^n
			break
		end
		
		RHS{2} = RHS_cell{j};
		if ~commutes(RHS)
			continue
		end
		
		[~, D{1}, D{2}, D{3}] = simdiag(RHS{1}, RHS{2}, LHS);
		
		for i = 1:runs+1
			d{i} = real(diag(D{i}));
			assert( all( imag(diag(D{i})) < 10e-5 ) );
		end

		for i = 1:runs
			const(i) = max(d{runs+1}-d{i});
			d{i} = d{i} + const(i);
		end
		
		rhs = d{1};
		for i = 2:runs
			rhs = min(rhs,d{i});
		end
		
		if all( (rhs - d{runs+1}) < 10e-5 )
			fprintf("LHS = ");
			print_envelope(coeffs(k,:), coeffs(j,:), int8(const), terms);
			count = count + 1;
		end
		%{
		if (mod(k,10) == 0) && (j == k+1)
			fprintf("%d\n",k);
		end
		%}
	end
end

fprintf('Found %d envelopes!\n', count);

function flag = commutes(RHS)
	flag = all( abs( RHS{1} * RHS{2} - RHS{2} * RHS{1} ) < 10e-5, 'all' );
end

function [] = print_envelope(quad1, quad2, const, terms)
	fprintf('min(');  print_quad(quad1, const(1), terms);
	fprintf(',');     print_quad(quad2, const(2), terms);
	fprintf(' )\n');
end

function [] = print_quad(quad, const, terms)
	for i = 1:numel(quad)
		if quad(i) == 1
			fprintf(' + %s',terms{i});
		elseif quad(i) == -1
			fprintf(' - %s',terms{i});
		end
	end
	if const > 0
		fprintf(' + %d',const);
	elseif const < 0
		fprintf(' - %d',-const);
	end
end
