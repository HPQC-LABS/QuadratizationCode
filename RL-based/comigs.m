clear;
global    stop_condition coef found num_COMIGs

% hyperparameters

fileID = 'LHS.txt';
num_cases = 1;
new_input = true;
seed = 'shuffle';

param.c_conflict = 1;
param.c_distance = 1;
param.c_entropy  = 0.25;

param.n_sessions = 100;
param.t_max = 150;

param.c_num_of_terms = 1;
param.c_strength = 2;
param.c_accuracy = 10000;

percentile = 70;
max_comigs = 16;

% -----------------------------------------------------------------------------------------

% Run for multiple LHSs

best_comig_cost = 1000;

%LHS_string = 'b1b2b3 + b3b4b5 + b4b5b6/n0';
file1 = fopen(fileID,'r');
LHS_string = fscanf(file1,'%s');
%close(file1);

aux = LHS_string(end) - '0';
LHS_string = LHS_string(1:end-1);

for Case = 0:num_cases
    % Read new input
    
    if new_input
        
        [init_training, LHS, allbits, reset_state, n] = pretrain(aux, LHS_string);
		
        coeffs_size = n*(n+1)/2;
        n_actions = 2*coeffs_size + 1;   % number of actions
        
        LHS_size = 2^(n-aux);
        
        LHS = reshape(LHS,LHS_size,1);
        
        allbits = allbits';
        
        reset_state = repmat(reset_state',1,param.n_sessions);
        
        one_hot = zeros(coeffs_size,n_actions);
        for i = 1:coeffs_size
            one_hot(i,i) = 1;
            one_hot(i,coeffs_size + i) = -1;
		end
		
		cache.LHS = LHS;
		cache.allbits = allbits;
		cache.n = n;
		cache.aux = aux;
		cache.coeffs_size = coeffs_size;
		cache.reset_state = reset_state;
		cache.one_hot = one_hot;
		clear LHS allbits n aux coeffs_size reset_state one_hot
    end
    % create agent
    
    agent = patternnet(cache.coeffs_size);
    agent.trainparam.epochs = 1;
    agent.trainparam.showWindow = false;
    
    Size = size(init_training,1);
    train_data = init_training - cache.one_hot(:,1)';
    target = ones(Size,1);
    
    for action = 2:n_actions
        train_data = [train_data; init_training - cache.one_hot(:,action)' ];
        target = [target; ones(Size,1)*action ];
    end
    
    train_data = train_data';
    target = ind2vec(target',n_actions);
    
    %------------------------------------------------------------------------------------------------
    % Finding COMIGs
    
    best_num_of_comigs = 1000;
    for Run = 1:3
        coef = zeros(max_comigs,cache.coeffs_size);
        score = 0;
        
        for num_COMIGs = 1:max_comigs
            
            mask = true(LHS_size,1);
            mask_preserve = zeros(max_comigs,LHS_size);
            
            for i = 1:num_COMIGs-1
                RHS = reshape( rhs(reshape(coef(i,:),[],1), cache) , [] , 1);
                mask_preserve(i,:) = (RHS ~= cache.LHS)';
                mask = mask.*mask_preserve(i,:)';
            end
            mask = reshape(mask,[],1);
            mask = (mask == 1);
            if sum(mask) == 0
                
                break;
            end
            fprintf('states left = %d\n', sum(mask));
            
            best_run_score = -10000;
            stop_counter = 0;
            counter = 0;
            found = 0;
            while counter < 4
                stop_condition = 0;
                rng(seed);
                agent = train(agent,train_data,full(target)); % initialize agent
                for j = 0:100
                    %mean_reward_max_in = -10000
                    
                    [batch_states, batch_actions, batch_rewards, flag, best_run_score] = generate_session(agent, param, cache, mask, best_run_score);
                    
                    if flag
                        break;
                    end
                    
                    [elite_states, elite_actions] = select_elites(batch_states, batch_actions, batch_rewards, percentile, cache);
                    
                    elite_actions = ind2vec(elite_actions, n_actions); % Transform to hot-one
                    agent = train(agent,elite_states, full(elite_actions));
                    
                    %{
                    flag = False
                    if j % 10 == 0:
                        flag = True
                    %}
                    mean_reward = show_progress(batch_rewards, percentile);%, flag)
                    %{
                    if mean_reward > mean_reward_max_in:
                        mean_reward_max_in = mean_reward
                
                if (mean_reward_max_in / (c_distance + c_entropy) / t_max) > mean_reward_max:
                    mean_reward_max = mean_reward_max_in
                    i_max = [t_max, c_distance, c_entropy]
                    print(i_max)
                        
                    if np.mean(batch_rewards)> 50:
                        print("You Win!!! You may stop training now")
                %}
                end
                if (counter == 3) && (found == 0) && (stop_counter < 6)
                    counter = 2;
                    stop_counter = stop_counter + 1;
                end
                counter = counter + 1;
            end
            fprintf('Case: %d, Run: %d, Comig: %d\n', Case, Run, num_COMIGs);
            
            coef_temp = coef(num_COMIGs,:);
            [conflict_current, ~,~] = accuracy(reshape(coef_temp,[],1), cache, mask);
            coeffs_strength = sum(abs(coef_temp));
            coeffs_number = sum(coef_temp ~= 0);
            best_score = param.c_accuracy*(1 - conflict_current) - param.c_strength*coeffs_strength - param.c_num_of_terms*coeffs_number;
            score = score + best_score;
        end
        
        fprintf('num of comigs = %d\n', num_COMIGs);
        if num_COMIGs < best_num_of_comigs
            best_num_of_comigs = num_COMIGs;
            Best_score = score;
            best_coef = coef(1:num_COMIGs);
            %savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
            %LATER!
        elseif num_COMIGs == best_num_of_comigs
            if score > Best_score
                best_coef = coef(1:num_COMIGs);
                %savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
                %LATER!
            end
        end
    end
    %{
	comig_cost = best_num_of_comigs*2^aux;
    
    if comig_cost < best_comig_cost
        best_comig_cost = comig_cost;
        best_aux = aux;
    end
    %}
    fprintf('Done!');
    fprintf('best num of comigs = %d\n',best_num_of_comigs);
    
    % verify coef and get const_terms
    [verified, const_terms] = verify(cache, best_coef);
    
    %file1 = fopen(fileID,'a');
    
    %fprintf(file1,'\n\n%d-aux\n', aux);
    %{
    if any(min(rhs(best_coef'),[],2) ~= LHS')
        file1.write('Not verified!\n');
    else
        file1.write('\nVerified!\n');
    end
    
    for m = 1:size(best_coef,1)
        vec = best_coef(m,:);
        k = 0;
        for i = 1:n
            for j = i+1:n+1
                if vec(k) ~= 0
                    if vec(k) == 1
                        file1.write(' + b_{%d}b_{%d}' , i, j );
                    elseif vec(k) == -1
                        file1.write(' - b_{%d}b_{%d}' , i, j );
                    else
                        file1.write(' %+d b_{%d}b_{%d}' , vec(k), i, j );
                    end
                k = k + 1;
                end
            end
        end
    
        for i = 1:n+1
            if vec(k) ~= 0
                if vec(k) == 1
                    file1.write(' + b_{%d}' , i);
                elseif vec(k) == -1
                    file1.write(' - b_{%d}' , i);
                else
                    file1.write(' %+d b_{%d}' , vec(k), i );
                end
            end
            k = k + 1;
        end
        
        if floor(const_terms(m)) ~= 0
            file1.write(' %+d' , floor(const_terms(m)) );
        end
        file1.write('\n');
    
    file1.close();
    end
    file1 = open(fileID,"a");

    file1.write('\nbest aux: %d\n' , best_aux);
    file1.write('best comig cost: %d' , best_comig_cost);

    file1.close();
    %}
end


%-----------------------------------------------------------------------------------------------

% functions

function r = reward(conflict, dist, entropy, param)
    r_conflict = 1 - conflict.^(1/3);         % "preserving states" reward
    r_distance = 1 - dist.^(1/3);             % "staying close"     reward
    r_entropy  = entropy;                     % "being an explorer" reward
    
    r = param.c_conflict * r_conflict + param.c_distance * r_distance + param.c_entropy * r_entropy;
end

function RHS = rhs(s, cache)
    RHS = cache.allbits*s;
    RHS = RHS - min(RHS) + min(cache.LHS);
    
    for i=1:cache.aux
        RHS = min(RHS(1:2:end,:),RHS(2:2:end,:));
    end
end

function [conflict, distance, flag] = accuracy(s, cache, mask)
    RHS = rhs(s, cache);
    difference = RHS - cache.LHS;
    flag = (sum(difference < 0) ~= 0);
    
    conflict = sum(difference(mask,:) ~= 0) / sum(mask);
    
    difference = difference(difference < 0);
    difference = abs(difference);
    distance = sum(difference).*flag;
    distance = min(1, distance/sum(mask)); %LHS_size? or maybe changing at runtime?
end

function [states, actions, total_reward, flag, best_score] = generate_session(agent, param, cache, mask, best_score)
    global  stop_condition coef found  num_COMIGs
    states  = zeros(size(cache.reset_state,1),size(cache.reset_state,2),param.t_max);
    actions = zeros(param.n_sessions,param.t_max);
    total_reward = 0;
    mini = 100;
    %r_best = -1000 * np.ones((n_sessions,))
    s = cache.reset_state;
    %conflict,dist = accuracy(s)
    %old_reward = reward(conflict, dist, 0)
    
    for t = 1:param.t_max
		
		probs = agent(s);
        
        % choose actions w.r.t. probs
        c = cumsum(probs);
        u = rand(1,size(c,2));
        [~,a] = max(u < c);
        
        new_s = s + cache.one_hot(:,a);
		
        [conflict, dist, flag] = accuracy(new_s, cache, mask);
        
        entropy = -sum( probs.*log(probs)/log(cache.coeffs_size));
        
        %new_reward = reward(conflict, dist, entropy)
        %r = new_reward - old_reward
        %old_reward = new_reward
        
        r = reward(conflict, dist, entropy, param);
        r(flag) = r(flag) - 10;
        
        condition = ( flag == 0 );
        
        if any(condition)
            coeffs_strength = sum(abs(new_s(:,condition)),1);
            coeffs_number = sum( new_s(:,condition) ~= 0 , 1 );
            score = param.c_accuracy*(1 - conflict(condition)) - param.c_strength*coeffs_strength - param.c_num_of_terms*coeffs_number;
            
			[score_max, score_idx] = max(score);
			
			if score_max > best_score
				best_score = score_max;
				coeffs = new_s(:,condition)';
				coef(num_COMIGs,:) = coeffs(score_idx,:);
				stop_condition = 0;
				found = 1;
			end
		end
		
        if min(conflict) < mini
            mini = min(conflict);
        end
        
        states(:,:,t) = s;
        actions(:,t)  = a;
        total_reward = total_reward + r;
        
        s = new_s;
    end
    
    stop_condition = stop_condition + 1;
    flag = false;
    if stop_condition == 40
        flag = true;
    end
	%print(1-mini)
end
        
function [elite_states, elite_actions] = select_elites(states_batch,actions_batch,rewards_batch,percentile, cache)
    %{
    Select states and actions from games that have rewards >= percentile
    :cache states_batch: list of lists of states, states_batch[session_i][t]
    :cache actions_batch: list of lists of actions, actions_batch[session_i][t]
    :cache rewards_batch: list of rewards, rewards_batch[session_i][t]
    
    :returns: elite_states,elite_actions, both 1D lists of states and respective actions from elite sessions
    %}
    
    reward_threshold = prctile(rewards_batch, percentile);
    mask = (rewards_batch > reward_threshold);
    
    elite_states = reshape(states_batch(:,mask,:), cache.coeffs_size, []);
    elite_actions = reshape(actions_batch(mask,:), 1, []);
end

function mean_reward = show_progress(batch_rewards, percentile)
    
    mean_reward = mean(batch_rewards);
    threshold = prctile(batch_rewards, percentile);
    %print("mean reward = %.3f, threshold=%.3f"%(mean_reward, threshold))
    %Log.append([int(mean_reward), int(threshold)])
    
    fprintf("mean reward = %.3f, threshold = %.3f\n" , mean_reward, threshold);
    %{
    plt.figure(figsize=[8,4])
    plt.subplot(1,2,1)
    plt.plot(list(zip(*log))[0], label='Mean rewards')
    plt.plot(list(zip(*log))[1], label='Reward thresholds')
    plt.legend()
    plt.grid()
    
    plt.subplot(1,2,2)
    plt.hist(batch_rewards, range=reward_range);
    plt.vlines([np.percentile(batch_rewards, percentile)], [0], [100], label="percentile", color='red')
    plt.legend()
    plt.grid()

    plt.show()
    %}
end

