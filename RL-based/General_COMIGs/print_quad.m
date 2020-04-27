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
