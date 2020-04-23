function [] = print_envelope(quad1, quad2, const, terms)
	fprintf('min(');  print_quad(quad1, const(1), terms);
	fprintf(',');     print_quad(quad2, const(2), terms);
	fprintf(' )\n');
end