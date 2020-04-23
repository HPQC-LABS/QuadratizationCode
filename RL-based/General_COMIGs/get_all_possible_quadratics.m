function [allbits, terms_] = get_all_possible_quadratics(n)

	sigma = cell(4,1);
	sigma{1} = [0 1 ; 1 0];
	sigma{2} = [0 -1i ; 1i 0];
	sigma{3} = [1 0 ; 0 -1];
	sigma{4} = eye(2);

	if n==3
		terms = cell(1,1000);
		term_idx = 0;
		allbits = zeros(2^n,0);
		%All 3-qubit combinations of X,Z,I that are up to quadratic
		for i=1:4
			for j=1:4
				for k=1:4
					if (i==4||j==4||k==4) && (i+j+k~=12) %&& (i~=2)&&(j~=2)&&(k~=2)
						allbits = [allbits kron(sigma{i},kron(sigma{j},sigma{k}))]; %[x1x2,x1z2,x1x3,x1z3,...,x3,z3,1];
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
		end

		terms_ = cell(1,term_idx);
		for i = 1:term_idx
			terms_{i} = terms{i};
		end
	elseif n==4
		terms = cell(1,1000);
		term_idx = 0;
		allbits = zeros(2^n,0);
		%All 4-qubit combinations of X,Z,I that are up to quadratic
		for i=1:4
			for j=1:4
				for k=1:4
					for m=1:4
						if ( ( ((i==4)&&(j==4)) || ((i==4)&&(k==4)) || ((i==4)&&(m==4)) || ((j==4)&&(k==4)) || ((j==4)&&(m==4)) || ((k==4)&&(m==4)) ) && (i~=2)&&(j~=2)&&(k~=2)&&(m~=2)&&(i+j+k+m~=16) )
							allbits = [allbits kron(sigma{i},kron(sigma{j},kron(sigma{k},sigma{m})))]; %[x1x2,x1z2,x1x3,x1z3,...,x3,z3];
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
		end

		terms_ = cell(1,size(terms,2));
		for i = 1:size(terms,2)
			terms_{i} = terms{i};
		end
	end
end
