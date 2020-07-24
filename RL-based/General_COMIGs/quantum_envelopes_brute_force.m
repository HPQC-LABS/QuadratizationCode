function quantum_envelopes_brute_force(H, c, X)
%	Inputs:
%		H - a Hamiltonian to get quadratized
%		c - coefficient constant (coeffs for the quadratic terms in [-c,c])
%		X - optional variable, discard quads that preserve < X states
%

warning('off');

[LHS, n] = get_LHS(H);
[allbits, terms] = get_all_possible_quadratics(n);
param = initialize_sim_diag(LHS);
lhs = param.d;	rhs = zeros(size(lhs));
J = param.J;
K = param.K;
cutoff = param.cutoff;

allbits_unfolded = reshape( allbits, 2^(2*n), []);
LHS_allbits_unfolded = reshape( LHS  * allbits, 2^(2*n), []);
allbits_LHS  = cell2mat( mat2cell( allbits' * LHS, ones(1,size(allbits_unfolded,2))*2^n, 2^n )' );
allbits_LHS_unfolded = reshape( allbits_LHS, 2^(2*n), []);
commutator_unfolded = LHS_allbits_unfolded - allbits_LHS_unfolded;
null_space = null(commutator_unfolded,'r');

print(H);
%print_nullspace(null_space, terms);
%null_space = null_space(:,[1,3,4]); % restrict the nullspace

if isempty(null_space), fprintf('Nothing commutes with the LHS! The Null Space is empty.\n'), return, end
allbits_unfolded = allbits_unfolded * null_space;
allbits_size = size(allbits_unfolded,2);

fprintf('%d Hamiltonians commute with LHS\n',size(null_space,2));
if size(null_space,2) > 15, return, end
if ~ishermitian(LHS), fprintf("Error! The LHS is not Hermitian!\n"); return, end

%{
switch size(null_space,2)
	case 1, c = 15;
	case 2, c = 15;
	case 3, c = 15;
	case 4, c = 10;
	case 5, c = 5;
	case 6, c = 3;
	case 7, c = 2;
	case 8, c = 2;
	otherwise, c = 1;
end
%}
base = size(-c:c,2);
coeffs_size = allbits_size;

restartId = 0;
perCheck = 100000;

data_percentage = 1;
data_size = floor(data_percentage*base^coeffs_size);
progress_const = 100 * perCheck / data_size;

coeffs_all = []; %zeros(data_size,coeffs_size);
preserved  = zeros(1,100000);

t = cputime;
t_init = t;
for checkpoint = restartId : floor( (data_size-1)/perCheck )
	k = int64(checkpoint*perCheck) : min(data_size-1, int64((checkpoint+1)*perCheck) - 1 );
	coeffs = ndec2base(k,base,coeffs_size) - '0';
	coeffs(coeffs > 9) = coeffs(coeffs > 9) - 7; % this is a fix in case base > 10
	% Putting them in order 0, 1, -1, 2, -2, ...
	for i = 2:2*c
		if mod(i,2)
			coeffs(coeffs == i) = (coeffs(coeffs == i) + 1)/2;
		else
			coeffs(coeffs == i) = -coeffs(coeffs == i)/2;
		end
	end
	RHS_unfolded = allbits_unfolded * coeffs';
	RHS = reshape( RHS_unfolded , 2^n, []);

	coeffs_all(k + 1,:) = coeffs;

	for i = 1:size(coeffs,1)
		RHS_ = RHS(:, 2^n*(i-1)+1:2^n*i );
		
		% handle degenerate eigenspaces
		A = param.Q'*RHS_*param.Q;
		for it = 1:numel(J)
			j = J(it); k = K(it);
			rhs(j:k) = real(eig(A(j:k,j:k)));
		end
		
		rhs = rhs + max(lhs - rhs);
		mask_preserve = (abs(rhs - lhs) < 10e-5);
		preserved(checkpoint * perCheck + i) = sum(mask_preserve);
	end
	
	if mod(checkpoint,1) == 0
		fprintf('progress %.2f%%, restart id = %d, step time = %.2f, total time = %.0f, found(half) = %d, found(1/3) = %d, found(1/4) = %d\n',...
			min((checkpoint+1)*progress_const, 100), checkpoint, cputime - t,cputime - t_init, sum(preserved >= 2^n/2), sum(preserved >= 2^n/3), sum(preserved >= 2^n/4));
		t = cputime;
	end
end

for runs = 2:4
	if quantum_envelopes_matching(runs, n, LHS, coeffs_all, preserved, allbits_unfolded, terms, null_space, H, c, X)
		break;
	end
end
end

function [LHS, n] = get_LHS(H, alpha)
	n = max(strlength(H));
	if nargin == 1, alpha = ones(1,size(H,2)); end
	
	LHS = zeros(2^n);
	for t = 1:size(H,2)
		LHS = LHS + alpha(t) * get_term(H{t}, n);
	end
end

function print(H, alpha)
	if nargin == 1, alpha = ones(1,size(H,2)); end
	
	for i = 1:size(H,2)
		if alpha(i) == 1
			if i ~= 1
				fprintf(' + %s', H{i});
			else
				fprintf('%s', H{i});
			end
		else
			fprintf(' %2d%s', alpha(i), H{i});
		end
	end
	fprintf(':\n');
end

