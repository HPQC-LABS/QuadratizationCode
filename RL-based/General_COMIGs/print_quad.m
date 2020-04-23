function [] = print_quad(quad, const, terms)
	for i = 1:numel(quad)
		if quad(i) == 1
			fprintf(' + %s',terms{i});
		elseif quad(i) == -1
			fprintf(' - %s',terms{i});
		elseif quad(i) > 0
			fprintf(' + %d%s',quad(i),terms{i});
		elseif quad(i) < 0
			fprintf(' - %d%s',-quad(i),terms{i});
		end
	end
	if const > 0
		fprintf(' + %d',const);
	elseif const < 0
		fprintf(' - %d',-const);
	end
end