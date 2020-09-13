function [allbits, terms] = get_all_possible_quadratics(n, skip_pauli)
	if nargin < 2
		skip_pauli = 'n';
		file = ['all_possible_quadratics_', int2str(n), '.mat'];
		if exist(file, 'file')
			load(file, 'allbits', 'terms');
			return;
		end
	elseif (skip_pauli ~= 'x') && (skip_pauli ~= 'y') && (skip_pauli ~= 'z')
		fprintf(2,'Error! Wrong input on get_all_possible_quadratics(n, skip_pauli)\n');
	end
	
	sigma = cell(4,1);
	sigma{1} = [0 1 ; 1 0];
	sigma{2} = [0 -1i ; 1i 0];
	sigma{3} = [1 0 ; 0 -1];
	sigma{4} = eye(2);
	
	terms = cell(1,4^n);
	term_idx = 0;
	allbits = zeros(2^n,0);
	
	if n == 3
		for i=1:4
			for j=1:4
				for k=1:4
					allbits = [allbits kron(sigma{i},kron(sigma{j},sigma{k}))];
					term = char([]);
					if i~= 4
						term = [term, char(i+119), '1'];
					end
					if j~= 4
						term = [term, char(j+119), '2'];
					end
					if k~= 4
						term = [term, char(k+119), '3'];
					end
					term_idx = term_idx + 1;
					terms{term_idx} = term;
				end
			end
		end
	elseif n == 4
		for i=1:4
			for j=1:4
				for k=1:4
					for m=1:4
						allbits = [allbits kron(sigma{i},kron(sigma{j},kron(sigma{k},sigma{m})))];
						term = char([]);
						if i~= 4
							term = [term, char(i+119), '1'];
						end
						if j~= 4
							term = [term, char(j+119), '2'];
						end
						if k~= 4
							term = [term, char(k+119), '3'];
						end
						if m~= 4
							term = [term, char(m+119), '4'];
						end
						term_idx = term_idx + 1;
						terms{term_idx} = term;
					end
				end
			end
		end
	elseif n == 5
		for i=1:4
			for j=1:4
				for k=1:4
					for m=1:4
						for p=1:4
							allbits = [allbits kron(sigma{i},kron(sigma{j},kron(sigma{k},kron(sigma{m},sigma{p}))))];
							term = char([]);
							if i~= 4
								term = [term, char(i+119), '1'];
							end
							if j~= 4
								term = [term, char(j+119), '2'];
							end
							if k~= 4
								term = [term, char(k+119), '3'];
							end
							if m~= 4
								term = [term, char(m+119), '4'];
							end
							if p~= 4
								term = [term, char(p+119), '5'];
							end
							term_idx = term_idx + 1;
							terms{term_idx} = term;
						end
					end
				end
			end
		end
	end
	[allbits, terms] = KeepQuadratics(allbits, terms, skip_pauli);
	file = ['all_possible_quadratics_', int2str(n), '.mat'];
	save(file, 'allbits', 'terms');
end

function [allbits_, terms_] = KeepQuadratics(allbits, terms, skip_pauli)
	width = size(allbits,1);
	
	terms_ = cell(1,2); % Size unknown, it will fill up dynamically
	allbits_  = [];
	count = 0;
	for k = 1:size(terms,2)
		switch size(terms{k},2)
			case 2,	flag = (terms{k}(1) ~= skip_pauli);
			case 4, flag = (terms{k}(1) ~= skip_pauli) && (terms{k}(3) ~= skip_pauli);
			otherwise, flag = false;
		end
		if flag
			count = count + 1;
			terms_{count} = terms{k};
			allbits_ = [allbits_ allbits( :, (k-1)*width+1 : k*width ) ];
		end
	end
end

