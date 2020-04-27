function [count] = get_and_print_envelope(LHS, RHS, runs, Q, IE, count, coeffs, preserved, terms, null_space, n, file)
	
	[flag, const] = get_envelope(LHS, RHS, runs, Q, IE);
	if flag
		print_envelope(coeffs, preserved, int8(const), terms, null_space, n, file);
		print_envelope(coeffs, preserved, int8(const), terms, null_space, n, 1); %std output
		count = count + 1;
	end
end
