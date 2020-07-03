function [count] = get_and_print_envelope(LHS, RHS, runs, Q, IE, count, coeffs, preserved, terms, null_space, n, file)
	
	[flag, const] = get_envelope(LHS, RHS, runs, Q, IE);
	if flag
		print_envelope(coeffs, preserved, int8(const), terms, null_space, n, file);
		print_envelope(coeffs, preserved, int8(const), terms, null_space, n, 1); %std output
		count = count + 1;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [flag, const] = get_envelope(LHS, RHS, runs, Q, IE)
	D = cell(1,runs+1);
	d = cell(1,runs+1);
	const = zeros(1,runs);
	
	[~, D{:}] = simdiag(Q, IE, LHS, RHS{:});
	for i = 1:runs+1
		d{i} = real(diag(D{i}));
		assert( all( imag(diag(D{i})) < 10e-5 ) );
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = print_envelope(quads, pres, const, terms, null_space, n, file)
	fprintf(file,"LHS = min(");
	for i = 1:size(quads,1)
		print_quad(quads(i,:), const(i), terms, null_space, file);
		if i ~= size(quads,1)
			fprintf(file,',');
		else
			fprintf(file,' )');
		end
	end
	for i = 1:size(quads,1)
		fprintf(file, ', %d/%d', pres(i), 2^n);
	end
	fprintf(file,'\n');
end

function [] = print_quad(quad, const, terms, null_space, file)
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
