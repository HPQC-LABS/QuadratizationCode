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
