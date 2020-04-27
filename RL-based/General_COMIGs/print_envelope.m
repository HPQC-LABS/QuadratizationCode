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
