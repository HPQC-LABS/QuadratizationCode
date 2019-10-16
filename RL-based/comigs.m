clear;
global c_conflict c_distance c_entropy worst_score stop_condition coef found reset_state agent allbits LHS aux one_hot mask coeffs_size n_sessions c_accuracy c_strength  c_num_of_terms num_COMIGs

% hyperparameters

fileID = 'LHS.txt';
num_cases = 1;
new_input = true;
seed = 'shuffle';

c_conflict = 1;
c_distance = 1;
c_entropy  = 0.25;

n_sessions = 100;
t_max = 150;

c_num_of_terms = 1;
c_strength = 2;
c_accuracy = 10000;

percentile = 70;
max_comigs = 16;

% -----------------------------------------------------------------------------------------

% Run for multiple LHSs

best_comig_cost = 1000;

file1 = fopen(fileID,'r');
LHS_string = fscanf(file1,'%s');
fclose(file1);

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
        
        reset_state = repmat(reset_state',1,n_sessions);
        
        one_hot = zeros(coeffs_size,n_actions);
        for i = 1:coeffs_size
            one_hot(i,i) = 1;
            one_hot(i,coeffs_size + i) = -1;
        end
    end
    % create agent
    
    agent = patternnet([coeffs_size,2*coeffs_size]);
    agent.trainParam.epochs = 1;
    agent.trainParam.showWindow = false;
    
    Size = size(init_training,1);
    train_data = init_training - one_hot(:,1)';
    target = ones(Size,1);
    
    for action = 2:n_actions
        train_data = [train_data; init_training - one_hot(:,action)' ];
        target = [target; ones(Size,1)*action ];
    end
    
    train_data = train_data';
    target = ind2vec(target',n_actions);
    
    %------------------------------------------------------------------------------------------------
    % Finding COMIGs
    
    best_num_of_comigs = 1000;
    for Run = 1:3
        coef = zeros(max_comigs,coeffs_size);
        score = 0;
        
        for num_COMIGs = 1:max_comigs
            
            mask = true(LHS_size,1);
            mask_preserve = zeros(max_comigs,LHS_size,1);
            
            for i = 1:num_COMIGs-1
                RHS = reshape( rhs(reshape(coef(i),[],1)) , [] , 1);
                size(RHS)
                mask_preserve(i) = (RHS ~= LHS);
                mask = mask.*mask_preserve(i);
            end
            mask = reshape(mask,[],1);
            mask = (mask == 1);
            if sum(mask) == 0
                
                break;
            end
            fprintf('states left = %d\n', sum(mask));
            
            worst_score = -10000;
            stop_counter = 0;
            counter = 0;
            found = 0;
            while counter < 4
                stop_condition = 0;
                rng(seed);
                agent = train(agent,train_data,full(target)); % initialize agent
                for j = 0:100
                    %mean_reward_max_in = -10000
                    
                    [batch_states,batch_actions,batch_rewards, flag] = generate_session(t_max);
                    
                    if flag
                        break;
                    end
                    
                    [elite_states, elite_actions] = select_elites(batch_states, batch_actions, batch_rewards, percentile);
                    
                    elite_actions = ind2vec(elite_actions, n_actions);
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
            [conflict_current, ~,~] = accuracy(reshape(coef_temp,[],1));
            coeffs_strength = sum(abs(coef_temp));
            coeffs_number = sum(coef_temp ~= 0);
            best_score = c_accuracy*(1 - conflict_current) - c_strength*coeffs_strength - c_num_of_terms*coeffs_number;
            score = score + best_score;
        end
        
        fprintf('num of comigs = %d\n', num_COMIGs);
        if num_COMIGs < best_num_of_comigs
            best_num_of_comigs = num_COMIGs;
            Best_score = score;
            best_coef = coef(1:num_COMIGs);
            %savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
            %LATER!!!!!y
        elseif num_COMIGs == best_num_of_comigs
            if score > Best_score
                best_coef = coef(1:num_COMIGs);
                %savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
                %LATER!!!!
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
    [verified, const_terms] = verify(n, aux, LHS, allbits, best_coef);
    
    file1 = fopen(fileID,'a');
    
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

function state = step(state, action)
	global one_hot
	state = state + one_hot(:,action);
end

function r = reward(conflict, dist, entropy)
    global c_conflict c_distance c_entropy
    
    r_conflict = 1 - conflict.^(1/3);            % "preserving states" reward
    r_distance = 1 - dist.^(1/3);                % "staying close"     reward
    r_entropy  = entropy;                        % "being an explorer" reward
    
    r = c_conflict * r_conflict + c_distance * r_distance + c_entropy * r_entropy;
end

function RHS = rhs(s)
    global allbits LHS aux
    RHS = allbits*s;
    RHS = RHS - min(RHS) + min(LHS);
    
    for i=1:aux
        RHS = min(RHS(1:2:end,:),RHS(2:2:end,:));
    end
end

function [conflict, distance, flag] = accuracy(s)
    global LHS mask
    RHS = rhs(s);
    difference = RHS - LHS;
    flag = (sum(difference < 0) ~= 0);
    
    conflict = sum(difference(mask,:) ~= 0) / sum(mask);
    
    difference = difference(difference < 0);
    difference = abs(difference);
    distance = sum(difference).*flag;
    distance = min(1, distance/sum(mask)); %LHS_size? or maybe changing at runtime?
end

function [states, actions, total_reward, flag] = generate_session(t_max)
    global worst_score stop_condition coef found reset_state agent coeffs_size n_sessions c_accuracy c_strength  c_num_of_terms num_COMIGs
    states  = zeros(size(reset_state,1),size(reset_state,2),t_max);
    actions = zeros(n_sessions,t_max);
    total_reward = 0;
    mini = 100;
    %r_best = -1000 * np.ones((n_sessions,))
    s = reset_state;
    %conflict,dist = accuracy(s)
    %old_reward = reward(conflict, dist, 0)
    
    for t = 1:t_max
        
        probs = agent(s);
        
        % choose actions w.r.t. probs
        c = cumsum(probs);
        u = rand(1,size(c,2));
        [~,a] = max(u < c);
        
        new_s = step(s,a);
        
        [conflict, dist, flag] = accuracy(new_s);
        
        entropy = -sum( probs.*log(probs)/log(coeffs_size));
        
        %new_reward = reward(conflict, dist, entropy)
        %r = new_reward - old_reward
        %old_reward = new_reward
        
        r = reward(conflict, dist, entropy);
        r(flag) = r(flag) - 10;
        
        condition = ( flag == 0 );
        
        if any(condition)
            coeffs_strength = sum(abs(new_s(:,condition)),1);
            coeffs_number = sum( new_s(:,condition) ~= 0 , 1 );
            score = c_accuracy*(1 - conflict(condition)) - c_strength*coeffs_strength - c_num_of_terms*coeffs_number;
            score_condition = ( score > worst_score );
            
            if any(score_condition)
                stop_condition = 0;
                found = 1;
                %flag = True
                %for k in range(len(good)):
                %    if (good[k] == new_s[:,temp]).all():
                %        flag = False
                %if flag:
                %   good.append(new_s[:,temp])
                %print(rhs(new_s[:,temp]))
                coeffs = new_s(:,condition)';
                coeffs = coeffs(score_condition,:);
                score = score(score_condition);
                %print(j,max(score),worst_score)
                for i = 1:size(coeffs,1)
                    if score(i) > worst_score
                        coef(num_COMIGs,:) = coeffs(i,:);
                        
                        coef_temp = coeffs(i,:);
                        [conflict_current, ~, ~] = accuracy(coef_temp');
                        coeffs_strength = sum(abs(coef_temp));
                        coeffs_number = sum(coef_temp ~= 0);
                        
                        score_current = c_accuracy*(1 - conflict_current) - c_strength*coeffs_strength - c_num_of_terms*coeffs_number;
                        if score_current < worst_score
                            worst_score = score_current;
                        end
                    end
                end
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
        
function [elite_states, elite_actions] = select_elites(states_batch,actions_batch,rewards_batch,percentile)
    global coeffs_size
    %{
    Select states and actions from games that have rewards >= percentile
    :param states_batch: list of lists of states, states_batch[session_i][t]
    :param actions_batch: list of lists of actions, actions_batch[session_i][t]
    :param rewards_batch: list of rewards, rewards_batch[session_i][t]
    
    :returns: elite_states,elite_actions, both 1D lists of states and respective actions from elite sessions
    %}
    
    reward_threshold = prctile(rewards_batch, percentile);
    mask = (rewards_batch > reward_threshold);
    
    elite_states = reshape(states_batch(:,mask,:), coeffs_size, []);
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

