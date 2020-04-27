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
