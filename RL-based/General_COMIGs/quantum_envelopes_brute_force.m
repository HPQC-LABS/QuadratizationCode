%function quantum_envelopes_brute_force(alpha)

c = 5;
H = {'XZZY', 'YYZY', 'XXZY', 'YXZY'};
alpha = [1, 1, 1, 1, 1];
N_of_terms = size(H,2);
n = max(strlength(H));

[allbits, terms] = get_all_possible_quadratics(n);

LHS = zeros(2^n);
for t = 1:N_of_terms
	LHS = LHS + alpha(t) * get_term(H{t}, n);
end

[Q,IE] = initialize_sim_diag(LHS);

%allbits = sparse(allbits);
allbits_unfolded = reshape( allbits, 2^(2*n), []);

LHS_allbits_unfolded = reshape( LHS  * allbits, 2^(2*n), []);
allbits_LHS  = cell2mat( mat2cell( allbits' * LHS, ones(1,size(allbits_unfolded,2))*2^n, 2^n )' );
allbits_LHS_unfolded = reshape( allbits_LHS, 2^(2*n), []);

commutator_unfolded = LHS_allbits_unfolded - allbits_LHS_unfolded;

null_space = null(commutator_unfolded,'r');
if isempty(null_space), fprintf('Nothing commutes with the LHS! The Null Space is empty.\n'), return, end
allbits_unfolded = allbits_unfolded * null_space;
allbits_size = size(allbits_unfolded,2);

fprintf('%d\n',size(null_space,2));
if ~ishermitian(LHS), fprintf("Error! The LHS is not Hermitian!\n"), return, end

%{
[U,V] = eig(LHS);

V(1,1) = V(1,1) + 2;
V(2,2) = V(2,2) + 2;
LHS_prime = U*V*U';

dist_threshold = 6;
%}

%{
    for checkpoint = restartId : floor( (((2*N1+1)*(2*N2+1))^N_of_terms(T)-1)/perCheck )
        for k = checkpoint*perCheck : min( ((2*N1+1)*(2*N2+1))^N_of_terms(T) - 1, (checkpoint+1)*perCheck - 1)
            alpha = zeros(N_of_terms(T));
            
            
            RHS = zeros(2^n,2^n);
            for t = 1:N_of_terms(T)
                RHS = RHS + alpha(t)*kron(sigma{H{T}{t}(1)-87},...
                    kron(sigma{H{T}{t}(2)-87},sigma{H{T}{t}(3)-87}));
            end
            
            if (max(max(abs(RHS))) < 1e-5) || ( ~ishermitian(RHS) ) % need some coeff non-zero
                    continue;
            end
            
            
            coeffsQ = ones(allbits_size,1);
            numbers = int64(1:2^n);
            
            for i = int64(1): int64(2^(2^n)-1)
                indexlist = bitget(i, 2^n : -1: 1);
                indexlist = indexlist.*numbers;
                indexlist = indexlist(indexlist~=0);
                
                A = allbits(indexlist,:)';
                A = reshape(A,2^n,[]);
                B = A(:,1:allbits_size);
                for m = 1:size(indexlist,2)-1
                    B = [B;A(:,1+m*allbits_size:(m+1)*allbits_size)];
                end
                
                coeffsQ = B\reshape(RHS(indexlist,:)',[],1);
                RHS = allbits*kron(coeffsQ,eye(2^n));
                if ~ishermitian(RHS)
                    continue;
                end
                
                [V, d] = eig(RHS);
                RHS_spectrum = uniquetol( diag(d) , 1e-5 );
                RHS_gs = V(:, diag(d)==min(diag(d)) );
                Delta_E = abs( LHS_spectrum(1) - RHS_spectrum(1) );
                
                if ( ( Delta_E < 1e-5 ) && ( max(coeffsQ) > 1e-5 ) )
                    r = rank( [LHS_gs, RHS_gs] );
                    hasQuad(k+1) = 1;
                    if ( size(LHS_gs,2) + size(RHS_gs,2) - r ) == 0
                        temp1 = [temp1, coeffsQ];
                    elseif  ( ( r == size(RHS_gs,2) ) && ( size(LHS_gs,2) <= r ) )
                        temp3 = [temp3, coeffsQ];
                    else
                        temp2 = [temp2, coeffsQ];
                    end
                    for m=2:min( size(LHS_spectrum,2),size(RHS_spectrum,2) )
                        if abs( LHS_spectrum(m) - RHS_spectrum(m) ) > 1e-5
                            break;
                        end
                    end
                    if m > 2
                        temp4 = [temp4, coeffsQ];
                    end
                end
            end
            quadratisations{k+1}{1} = alpha;
            %remove multiple entries
            if ~sum(sum(real(temp1)~=temp1))
                temp1 = uniquetol(temp1',1e-5,'ByRows',true)';
            end
            if ~sum(sum(real(temp2)~=temp2))
                temp2 = uniquetol(temp2',1e-5,'ByRows',true)';
            end
            if ~sum(sum(real(temp3)~=temp3))
                temp3 = uniquetol(temp3',1e-5,'ByRows',true)';
            end
            if ~sum(sum(real(temp4)~=temp4))
                temp4 = uniquetol(temp4',1e-5,'ByRows',true)';
            end
            quadratisations{k+1}{2} = temp1;
            quadratisations{k+1}{3} = temp2;
            quadratisations{k+1}{4} = temp3;
            quadratisations{k+1}{5} = temp4;
        end
        
        for k = checkpoint*perCheck : min( ((2*N1+1)*(2*N2+1))^N_of_terms(T) - 1, (checkpoint+1)*perCheck - 1)
            alpha = zeros(N_of_terms(T));
            for t = 1:N_of_terms(T)
                temp = mod(floor( k / ((2*N1+1)*(2*N2+1))^(t-1) ), (2*N1+1)*(2*N2+1) );
                a = floor(temp/(2*N2+1))-N1;
                b = mod(temp,2*N2+1) - N2;
                alpha(t) = a + 1i*b;
            end
            if norm(alpha) < 1e-5 % need some coeff non-zero
                continue;
            end
            
            if hasQuad(k+1) == 0
                fprintf(fileid{5}, "No quadratisation for (%d%+di)%c1%c2%c3", real(alpha(1)), imag(alpha(1)), H{T}{1}(1), H{T}{1}(2), H{T}{1}(3));
                for m = 2:N_of_terms(T)
                    fprintf(fileid{5}, " + (%d%+di)%c1%c2%c3", real(alpha(m)), imag(alpha(m)), H{T}{m}(1), H{T}{m}(2), H{T}{m}(3));
                end
                fprintf(fileid{5}, "\n");
            else
                for tempm=1:4
                    temp = quadratisations{k+1}{tempm+1};
                    if size(temp,2)
                        plus_flag = false;
                        for m = 1:N_of_terms(T)
                                if plus_flag
                                    fprintf(fileid{tempm}, " + ");
                                end
                                fprintf(fileid{tempm}, "(%d%+di)%c1%c2%c3", real(alpha(m)), imag(alpha(m)), H{T}{m}(1), H{T}{m}(2), H{T}{m}(3));
                                plus_flag = true;
                        end
                        fprintf(fileid{tempm}, " has quadratisations:\n");
                        %print possible quadratisations
                        for count=1:size(temp,2)
                            m = 0;
                            plus_flag = false;
                            for i=1:4
                                for j=1:4
                                    for l=1:4
                                        if (i==4||j==4||l==4)
                                            m = m+1;
                                            if norm(temp(m,count)) > 1e-5
                                                if imag(temp(m,count)) > 1e-5
                                                    if plus_flag
                                                        fprintf(fileid{tempm}, " + ");
                                                    end
                                                    fprintf(fileid{tempm}, "(%.1f%+.1fi)", real(temp(m,count)), imag(temp(m,count)));
                                                else
                                                    fprintf(fileid{tempm}, " %+ .1f", real(temp(m,count)));
                                                end
                                                plus_flag = true;
                                                if(i~=4)
                                                    fprintf(fileid{tempm}, "%c1", char(87+i));
                                                end
                                                if(j~=4)
                                                	fprintf(fileid{tempm}, "%c2", char(87+j));
                                                end
                                                if(l~=4)
                                                	fprintf(fileid{tempm}, "%c3", char(87+l));
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            fprintf(fileid{tempm}, "\n");
                        end
                        fprintf(fileid{tempm}, "\n");
                    end
                end
            end
        end
        fprintf('progress %.3f%%, restart id = %d\n', min((checkpoint+1)*perCheck/(((2*N1+1)*(2*N2+1))^N_of_terms(T) - 1), 1) * 100, checkpoint);
	end

toc
warning('on', 'all');

%}

coeffs_range = -c:c;
base = size(coeffs_range,2);
init = int2str(abs(coeffs_range(1)));
coeffs_size = allbits_size;

restartId = 0;
perCheck = 10000;

data_percentage = 1;
data_size = floor(data_percentage*base^coeffs_size);
progress_const = 100 * perCheck / data_size;

coeffs_all = []; %zeros(data_size,coeffs_size);
preserved  = zeros(1,100000);

t = cputime;
t_init = t;
for checkpoint = restartId : floor( (data_size-1)/perCheck )
	k = int64(checkpoint*perCheck) : min( floor(data_size-1), int64((checkpoint+1)*perCheck) - 1 );
	%k = randperm( data_size-1 , perCheck); % use this for random sample
	coeffs = ndec2base(k,base,coeffs_size) - '0';%init;
	coeffs(coeffs > 9) = coeffs(coeffs > 9) - 7; % this is a fix in case base > 10
	for i = 2:2*c
		if mod(i,2)
			coeffs(coeffs == i) = (coeffs(coeffs == i) + 1)/2;
		else
			coeffs(coeffs == i) = -coeffs(coeffs == i)/2;
		end
	end
	%coeffs = sparse(coeffs);
	%{
	%commutator_unfolded_ = commutator_unfolded * coeffs';
	
	%LHS_RHS = coeffs * LHS_allbits_unfolded';
	%RHS_LHS = coeffs * allbits_LHS_unfolded;
	
	%LHS_RHS = reshape( LHS_RHS, 2^n, []);
	%RHS_LHS = reshape( RHS_LHS, [], 2^n);
	
	%{
	RHS_unfolded = allbits_unfolded * coeffs';
	RHS = reshape( RHS_unfolded , 2^n, []);
	
	LHS_RHS = LHS  * RHS;   % LHS  * RHS: size 8x8     * 8x80000 = 8     * 80000
	RHS_LHS = RHS' * LHS;   % RHS' * LHS: size 80000x8 * 8x8     = 80000 * 8
	%}
	
	%LHS_RHS = reshape( permute( reshape(LHS_RHS,2^n,2^n,[]), [1,3,2] ), [], 2^n );
	
	%{
	flag = any( abs( LHS_RHS - RHS_LHS ) > 10e-5 , 2);
	for i = 1:n
		flag = flag(1:2:end) + flag(2:2:end);
	end
	%}
	
	%dist = cellfun(@(x) norm(LHS_prime - x),RHS_cell);
	%}
	
	RHS_unfolded = allbits_unfolded * coeffs';
	RHS = reshape( RHS_unfolded , 2^n, []);
	
	coeffs_all(k + 1,:) = coeffs;
	
	for i = 1:size(coeffs,1)
		RHS_ = RHS(:, 2^n*(i-1)+1:2^n*i );
		[~,D1,D2] = simdiag(Q, IE, LHS, RHS_);
		
		lhs = real(diag(D1));  assert( all( imag(diag(D1)) < 10e-5 ) );
		rhs = real(diag(D2));  assert( all( imag(diag(D2)) < 10e-5 ) );

		rhs = rhs + max(lhs - rhs);
		mask_preserve = (abs(rhs - lhs) < 10e-5);
		preserved(checkpoint * perCheck + i) = sum(mask_preserve);
	end

	%{
		fprintf('Nice! %d\n',size(coeffs_all,1));
		[~,D1,D2] = simdiag(LHS,allbits*kron(coeffs_all(end,:)',eye(8)));
		d1 = diag(D1);
		d2 = diag(D2);
		const = max(d1-d2);
		d2 = d2 + const;
		mask_preserve = (d2-d1 < 10e-5);
		%}
	%{
	
	if ~all(flag)
		idx = find(flag == 0);
		coeffs_new = coeffs(idx,:);
		coeffs_all = [coeffs_all ; coeffs_new];
		
		RHS_unfolded = allbits_unfolded * coeffs_new';
		RHS = reshape( RHS_unfolded , 2^n, []);
		
		for i = 1:numel(idx)
			RHS_ = RHS(:, 8*i-7:8*i );
			[~,D1,D2] = simdiag(LHS, RHS_);
			
			lhs = real(diag(D1));  assert( all( imag(diag(D1)) < 10e-5 ) );
			rhs = real(diag(D2));  assert( all( imag(diag(D2)) < 10e-5 ) );

			const = max(lhs - rhs);
			rhs = rhs + const;
			mask_preserve = (abs(rhs - lhs) < 10e-5);
			%preserved{checkpoint}(pre_idx) = sum(mask_preserve);
			preserved(pre_idx) = sum(mask_preserve);
			pre_idx = pre_idx + 1;
		end

		%{
		fprintf('Nice! %d\n',size(coeffs_all,1));
		[~,D1,D2] = simdiag(LHS,allbits*kron(coeffs_all(end,:)',eye(8)));
		d1 = diag(D1);
		d2 = diag(D2);
		const = max(d1-d2);
		d2 = d2 + const;
		mask_preserve = (d2-d1 < 10e-5);
		%}
	end
	%}
	%{
	if min(dist) < dist_threshold
		coeffs_all = [coeffs_all ; coeffs(dist < dist_threshold,:)];
		%fprintf('Nice! %d\n',size(coeffs_all,1));
	end

	conflicts_percent = mean( RHS ~= LHS , 2 ) * 100; % percentage of overall conflicts
	index_good = (conflicts_percent <= conflicts_threshold);

	difference = RHS - LHS;
	flag = (sum(difference < 0,2) ~= 0);

	index_good = logical(index_good.*flag);
	%input(checkpoint*perCheck+1:min( floor(data_size), int64((checkpoint+1)*perCheck)),:) = coeffs; %for sampling
	%target(checkpoint*perCheck+1:min( floor(data_size), int64((checkpoint+1)*perCheck))) = index_good; %for sampling
	if any(index_good)
		good_coeffs = [good_coeffs; coeffs(index_good,:)];
		good_prcntg = [good_prcntg; conflicts_percent(index_good)];
		good_count = good_count + sum(index_good);
	end
	%}
	if mod(checkpoint,1) == 0
		fprintf('progress %.2f%%, restart id = %d, step time = %.2f, total time = %.0f, found4+ = %d, found6+ = %d\n',...
			min((checkpoint+1)*progress_const, 100), checkpoint, cputime - t,cputime - t_init, sum(preserved >= 4), sum(preserved >= 6));
		t = cputime;
	end

	%{
	if (size(coeffs_all,1) > 1500) || (cputime - t_init > 60)
		break;
	end
	%}
end

%{
input = coeffs_all(randperm(size(coeffs_all,1),min(1000,size(coeffs_all,1))),:);

accuracy_percent = mean( RHS == LHS , 2 ) * 100; % percentage of overall conflicts
index_good = (accuracy_percent >= 60);
if sum(index_good) == 0
	index_good = (accuracy_percent >= 40);
end

temp = good_coeffs(index_good,:);
reset_state = temp(randperm(size(temp,1),1),:);
%save('data.mat','input','LHS','allbits','reset_state');

reset_state = input(1,:);
%}

for runs = 2:5
	if quantum_envelopes_matching(runs, n, LHS, coeffs_all, preserved, allbits_unfolded, terms, null_space, alpha)
		break;
	end
end

