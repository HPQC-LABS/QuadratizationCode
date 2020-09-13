global allbits_unfolded LHS perCheck coeffs_all preserved data_size base coeffs_size n

H = {'ZZX'};
alpha = [1, 1];
N_of_terms = size(H,2);
n = 3;

sigma = cell(4,1);
sigma{1} = [0 1 ; 1 0];
sigma{2} = [0 -1i ; 1i 0];
sigma{3} = [1 0 ; 0 -1];
sigma{4} = eye(2);

if n==3
	terms = cell(1,18);
	term_idx = 0;
	allbits = zeros(2^n,0);
	%All 3-qubit combinations of X,Z,I that are up to quadratic
	for i=1:4
		for j=1:4
			for k=1:4
				if ( (i==4||j==4||k==4)&&(i~=2)&&(j~=2)&&(k~=2)&&(i+j+k~=12) )
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
	allbits_size = size(allbits,2)/2^n;

	LHS = zeros(2^n,2^n);
	for t = 1:N_of_terms
		LHS = LHS + alpha(t)*kron(sigma{H{t}(1)-87},...
			kron(sigma{H{t}(2)-87},sigma{H{t}(3)-87}));
	end
	
elseif n==4
	allbits = zeros(2^n,0);
	%All 4-qubit combinations of X,Z,I that are up to quadratic
	for i=1:4
		for j=1:4
			for k=1:4
				for m=1:4
					if ( ( ((i==4)&&(j==4)) || ((i==4)&&(k==4)) || ((i==4)&&(m==4)) || ((j==4)&&(k==4)) || ((j==4)&&(m==4)) || ((k==4)&&(m==4)) ) && (i~=2)&&(j~=2)&&(k~=2)&&(m~=2)&&(i+j+k+m~=16) )
						allbits = [allbits kron(sigma{i},kron(sigma{j},kron(sigma{k},sigma{m})))]; %[x1x2,x1z2,x1x3,x1z3,...,x3,z3,1];
					end
				end
			end
		end
	end
	allbits_size = size(allbits,2)/2^n;

	LHS = zeros(2^n,2^n);
	for t = 1:N_of_terms
		LHS = LHS + alpha(t)*kron(sigma{H{t}(1)-87},...
			kron(sigma{H{t}(2)-87},kron(sigma{H{t}(3)-87},sigma{H{t}(4)-87})));
	end
end

allbits = sparse(allbits);
allbits_unfolded = reshape( allbits, 2^(2*n), []);

if ~ishermitian(LHS)
	fprintf("ERROR!!!!!! LHS is not Hermitian!");
end

coeffs_range = -1:1;
coeffs_size = allbits_size;
base = size(coeffs_range,2);
init = int2str( (base-1)/2 );

restartId = 0;
perCheck = 10000;

data_percentage = 1;
conflicts_threshold = 60;

coeffs_all = cell(1,100);
preserved = cell(1,100);

data_size = floor(data_percentage*base^coeffs_size);
progress_const = 100 * perCheck / data_size;
t = cputime;
t_init = t;

numCores = nproc();

%checkpoint = restartId : floor( (data_size-1)/perCheck )

pararrayfun(numCores, @find_envelopes, 0:99);

%{
fprintf('progress %.4f%%, restart id = %d, step time = %.3f, total time = %.3f, found = %d\n',...
	min((checkpoint+1)*progress_const, 100), checkpoint, cputime - t,cputime - t_init, pre_idx);
t = cputime;
%}

