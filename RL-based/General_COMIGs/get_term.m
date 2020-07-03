function term = get_term(str, n)
	sigma = cell(4,1);
	sigma{1} = [0 1 ; 1 0];
	sigma{2} = [0 -1i ; 1i 0];
	sigma{3} = [1 0 ; 0 -1];
	sigma{4} = eye(2);
	
	if str(1) == 'I'
		term = sigma{4};
	else
		term = sigma{str(1) - 87};
	end
	for i = 2:n
		if str(i) == 'I'
			term = kron(term, sigma{4});
		else
			term = kron(term, sigma{str(i) - 87});
		end
	end
end
