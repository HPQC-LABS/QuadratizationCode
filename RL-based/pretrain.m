function [input,LHS,allbits,reset_state,n] = pretrain(aux,LHS_string)
	
    b_ind = regexp(LHS_string,'b*');
    numbers = LHS_string(b_ind+1) - '0';
    n = max(numbers) + aux;
	if max(numbers) == 9
		maxnum = 0;
		b_ind = regexp(LHS_string,'b1');
		for i = 1:size(b_ind,2)
			num = sscanf(LHS_string(b_ind(i)+1:end),'%d');
			if num > maxnum
				maxnum = num;
			end
		end
		if maxnum > 9
			n = maxnum + aux;
		end
	end
	
    coeffs_range = -1:1;
    allCombos = dec2bin(0:2^n-1) -'0';

    b = cell(n,1);
    for i = 1:n
        b{i} = allCombos(:,i);
    end
    
    LHS = 0;
    term = 1;
    if LHS_string(1) == '+' || LHS_string(1) == '-'
    	term = 0;
    end
    i = 1;
    while i <= size(LHS_string,2)
        switch LHS_string(i)
            case 'b'
				temp = sscanf(LHS_string(i+1:end),'%d');
                term = term.*b{temp(1)};
                i = i + 1 + floor( log(temp)/log(10) );
			case 'x'
				temp = sscanf(LHS_string(i+3:end),'%d');
                term = term.*b{temp(1)};
                i = i + 4 + floor( log(temp)/log(10) );
            case '+'
                LHS = LHS + term;
                term = 1;
            case '-'
                LHS = LHS + term;
                term = -1;
            case '\'
                break;
            case num2cell('1':'9')
                temp = sscanf(LHS_string(i:end),'%d');
                term = term * temp;
                i = i + floor( log(temp)/log(10) ) ;
        end
        i = i + 1;
    end
    LHS = LHS + term;
    
    
    %LHS = b{1}.*b{2}.*b{3}.*b{4} + b{4}.*b{5}.*b{6}.*b{7};
    LHS = LHS';
    LHS = LHS(1:2^aux:end);
    
    allbits = [];
    for i = 1:n
        for j = i+1:n
            allbits = [allbits b{i}.*b{j}];
        end
    end

    for i = 1:n
        allbits = [allbits b{i}];
    end
	%allbits = [allbits ones(2^n,1)];
    allbits = allbits';

    coeffs_size = n*(n+1)/2;
    base = size(coeffs_range,2);
    init = int2str( (base-1)/2 );

    restartId = 0;
    perCheck = 1000;
    
    data_percentage = 1;
    conflicts_threshold = 60;
    
    data_size = floor(data_percentage*base^coeffs_size);
    progress_const = 100 * perCheck / data_size;
    good_coeffs = [];
    good_prcntg = [];
    good_count = 0;
    t = cputime;
    t_init = t;
    for checkpoint = restartId : floor( (data_size-1)/perCheck )
        k = int64(checkpoint*perCheck) : min( floor(data_size-1), int64((checkpoint+1)*perCheck) - 1 );
        %k = randperm( data_size-1 , perCheck); % use this for random sample
        coeffs = ndec2base(k,base,coeffs_size) - '0';%init;
        coeffs( coeffs == 2 ) = -1;
				
        RHS = coeffs*allbits;
        for i = 1:aux
            RHS = min(RHS(:,1:2:end),RHS(:,2:2:end));
        end
        RHS = RHS - min(RHS,[],2) + min(LHS);

        conflicts_percent = mean( RHS ~= LHS , 2 ) * 100; % percentage of overall conflicts
        index_good = (conflicts_percent <= conflicts_threshold);
        
        difference = RHS - LHS;
        flag = (sum(difference < 0, 2) ~= 0);
        
        index_good = logical(index_good.*flag);
        %input(checkpoint*perCheck+1:min( floor(data_size), int64((checkpoint+1)*perCheck)),:) = coeffs; %for sampling
        %target(checkpoint*perCheck+1:min( floor(data_size), int64((checkpoint+1)*perCheck))) = index_good; %for sampling
        if any(index_good)
            good_coeffs = [good_coeffs; coeffs(index_good,:)];
            good_prcntg = [good_prcntg; conflicts_percent(index_good)];
            good_count = good_count + sum(index_good);
        end
        %{
        if mod(checkpoint,100) == 0
            fprintf('progress %.4f%%, restart id = %d, step time = %.3f, total time = %.3f, good = %d\n',...
                min((checkpoint+1)*progress_const, 100), checkpoint,cputime - t,cputime - t_init,good_count);
            t = cputime;
        end
        %}
        if (size(good_coeffs,1) > 10000) || (cputime-t_init > 90)
            break;
        end
    end

    %fprintf('progress %.3f%%, restart id = %d, step time = %.3f, total time = %.3f, good = %d\n',...
    %    min((checkpoint+1)*progress_const, 100), checkpoint,cputime - t,cputime - t_init,good_count);
	
    input = good_coeffs(randperm(size(good_coeffs,1),min(1000,size(good_coeffs,1))),:);
	if ~isempty(good_coeffs)
		RHS = good_coeffs*allbits;
	    for i = 1:aux
		    RHS = min(RHS(:,1:2:end),RHS(:,2:2:end));
		end
		RHS = RHS - min(RHS,[],2) + min(LHS);
    
	    accuracy_percent = mean( RHS == LHS , 2 ) * 100; % 100 - percentage of overall conflicts
		index_good = (accuracy_percent >= 60);
		if sum(index_good) == 0
			index_good = (accuracy_percent >= (100 - conflicts_threshold));
		end

		temp = good_coeffs(index_good,:);
		reset_state = temp(randperm(size(temp,1),1),:);
		%save('data.mat','input','LHS','allbits','reset_state');
	end
end