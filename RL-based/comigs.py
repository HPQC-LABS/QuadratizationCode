import numpy as np
import re
from sklearn.neural_network import MLPClassifier
from sklearn.exceptions import ConvergenceWarning
import warnings
warnings.filterwarnings(action = 'ignore', category = ConvergenceWarning)
import time

# hyperparameters

profiling = 0   # 0: no prof, 1: function prof, 2: line prof
n_sessions = 100
t_max = 150
percentile = 70
max_runs = 20

filein  = "Quads/LHS.txt"
fileout = "Quads/quad"
#filein  = "LHS_monomial.txt"
#fileout = "quad_monomial"

#c_conflict = 1
#c_distance = 1
#c_entropy  = 0.25

#c_num_of_terms = 1
#c_strength = 1

############################# PRE-TRAINING ####################################


def dec2base(n, b = 2, padding = 0):
    num = np.base_repr(n, b)
    
    return num.zfill(padding)

def get_n(aux, LHS_string, latex):
    
    if latex:
        form = r"(b_{\d+)|(x_{\d+)"
        offset = 3
    else:
        form = r"(b\d+)|(x\d+)"
        offset = 1
    
    r1 = re.findall(form, LHS_string)
    r2 = [int(string[offset:]) for item in r1 for string in item if string != '']
    
    return max(r2) + aux

def get_b(n):
    form = "".join(['0',str(n),'b'])
    temp = [format(x,form) for x in range(2**n)]
    temp2 = [[int(item) for item in list(temp[it])] for it in range(len(temp))]
    allCombos = np.array(temp2)
    
    return np.transpose(allCombos)

def get_term(term_string, b, latex):
    
    if latex:
        form = r"(b_{\d+)|(x_{\d+)"
        offset = 3
    else:
        form = r"(b\d+)|(x\d+)"
        offset = 1
    
    if re.search(r"^\d+",term_string):
        temp = re.findall(r"^\d+",term_string)
        term = int(temp[0])
    else:
        term = 1
    
    bs = re.findall(form, term_string)
    b_idx = [int(string[offset:]) for item in bs for string in item if string != '']
    
    for idx in b_idx:
        term *= b[idx-1]
    
    return term

def get_LHS(b, aux, LHS_string, latex):
    
    if (LHS_string[0] != '-'):
        LHS_string = "".join(['+ ',LHS_string])
    
    LHS_split = LHS_string.split()
    
    LHS = 0
    for it in range(int(len(LHS_split)/2)):
        term = get_term(LHS_split[2*it + 1], b, latex)
        if (LHS_split[2*it] == '+'):
            LHS = LHS + term
        elif (LHS_split[2*it] == '-'):
            LHS = LHS - term
        else:
            print('Error reading LHS from text!')
    
    LHS = np.transpose(np.array(LHS))
    
    return LHS[::2**aux]

def pretrain(aux, LHS_string, pretrain_idx, latex = False):
    
    n = get_n(aux, LHS_string, latex)
    b = get_b(n)
    LHS = get_LHS(b, aux, LHS_string, latex)
    
    allbits = [ b[i]*b[j] for i in range(n) for j in range(i+1,n) ]
    
    for i in range(n):
        allbits.append(b[i])
    
    allbits.append(np.ones((2**n,))) # to include const terms in learning
    allbits = np.array(allbits)
    allbits = np.transpose(allbits)
    
    coeffs_size = int(n*(n+1)/2 + 1)
	
    base = 3    
    restartId = 0
    perCheck = 1000
    
    data_percentage = 1
    conflicts_threshold = 60
    
    data_size = np.floor(data_percentage*base**coeffs_size)
    vdec2base = np.vectorize(dec2base)
    found = 0
    good_count = 0
    t = time.time()
    t_init = t
    for checkpoint in range( restartId, int( (data_size-1)/perCheck ) ):
        k = np.arange(checkpoint*perCheck, min(int(data_size), (checkpoint+1)*perCheck) )
        
        temp = vdec2base(k, base, padding = coeffs_size)
        temp2 = [[int(item) for item in list(temp[it])] for it in range(len(temp))]
        coeffs = np.array(temp2)
        
        coeffs[ coeffs == 2 ] = -1
        
        RHS = np.transpose(rhs(np.transpose(coeffs), allbits, aux))

        conflicts_percent = np.mean( RHS != LHS , axis = 1) * 100 # percentage of overall conflicts
        index_good = conflicts_percent <= conflicts_threshold
        flag = (RHS < LHS).any(axis = 1)
        
        index_good = index_good * flag
        if index_good.any():
            if not found:
                good_coeffs = np.array(coeffs[index_good])
                good_prcntg = np.array(conflicts_percent[index_good])
                found = 1
            else:
                good_coeffs = np.append(good_coeffs, coeffs[index_good], axis = 0)
                good_prcntg = np.append(good_prcntg, conflicts_percent[index_good], axis = 0)
            good_count = good_prcntg.shape[0]
        
        if (good_count > 5000) or (time.time() - t_init > 10):
            break
    
    if good_coeffs.size:
		
        RHS = np.transpose(rhs(np.transpose(good_coeffs), allbits, aux))
        
        conflict_percent = np.mean( RHS != LHS , axis = 1 ) * 100
        index_good = (conflict_percent <= 40)
        if not index_good.any():
            index_good = (conflict_percent <= conflicts_threshold)
		
        temp = good_coeffs[index_good]
        if pretrain_idx == -1:
            pretrain_idx  = int(temp.shape[0] * np.random.rand())
        reset_state = temp[pretrain_idx]
    
    return good_coeffs, LHS, allbits, reset_state, n, pretrain_idx 


################################ CLASSES ######################################

class create_cache:
    
    def __init__(self, LHS, allbits, coeffs_size, n_actions, n, aux, reset_state, env):
        self.LHS_original = LHS.astype(int)
        self.LHS_ground = np.min(self.LHS_original)
        self.LHS = self.LHS_original - self.LHS_ground
        
        self.allbits = allbits.astype(np.float64)
        self.coeffs_size = coeffs_size
        self.log_n_actions = np.log(n_actions)
        self.n = n
        self.aux = aux
        self.reset_state = reset_state.astype(np.float64)
        self.env = env
    
    def add_mask(self, mask, mask_sum):
        self.mask = mask
        self.mask_sum = mask_sum

    def add_seed(self, seed):
        self.seed = seed

class environment:
    
    def __init__(self,one_hot):
        self.one_hot = one_hot.astype(np.float64)
    
    def reset(self, reset_state):
        self.state = reset_state
        return self.state
    
    def step(self, action):
        self.state = np.add(self.state, self.one_hot[:,action])
        return self.state


########################## AUXILIARY FUNCTIONS ################################


def reward(conflict, dist, entropy):
    r_conflict = 1 - conflict**(1/3)             # "preserving states" reward
    r_distance = 1 - dist**(1/3)                 # "staying close"     reward
    r_entropy  = entropy                         # "being an explorer" reward
    
    #return c_conflict * r_conflict + c_distance * r_distance + c_entropy * r_entropy - 10*flag
    return r_conflict + r_distance + 0.25 * r_entropy

def rhs(s, allbits, aux):
    RHS = allbits.dot(s)
    RHS.astype(int)
    
    for _ in range(aux):
        RHS  = np.minimum(RHS[::2],RHS[1::2])
    
    return RHS

def accuracy(s, cache, allbits, aux):
    
    RHS = rhs(s, allbits, aux)
    difference = RHS - cache.LHS
    
    conflict = np.sum( difference[cache.mask] != 0 , 0 ) / cache.mask_sum
    difference = np.where(difference < 0, difference, 0)
    distance = np.abs(np.sum(difference, 0))
    distance = np.minimum(1, distance/cache.mask_sum)
    
    flag = distance != 0
    
    return conflict, distance, flag

def generate_session(agent, cache, best_score, best_conflict, t_max = 50):
    global stop_condition, coef, run_number, seed_vec
    env = cache.env
    states,actions = [],[]
    total_reward = 0
    s = env.reset(cache.reset_state)
    allbits = cache.allbits
    aux = cache.aux
    
    for t in range(t_max):
        
        probs = agent.predict_proba(np.transpose(s))
        
        # choose actions w.r.t. probs
        c = probs.cumsum(axis = 1)
        u = np.random.rand(len(c),1)
        a = (u < c).argmax(axis = 1)
        
        new_s = env.step(a)
        
        conflict, dist, flag = accuracy(new_s, cache, allbits, aux)
        
        #entropy = enpy(np.transpose(probs), base = cache.n_actions)
        
        entropy = - np.sum( np.multiply(probs,np.log(probs)), 1) / cache.log_n_actions
        
        r = reward(conflict, dist, entropy)
        r[flag] -= 10
        
        if not all(flag):
            condition = np.logical_not(flag)
            conflict_new = conflict[condition]
            min_idx = np.argmin(conflict_new)
            if conflict_new[min_idx] < best_conflict:
                
                best_conflict = conflict_new[min_idx]
                coeff = new_s[:,condition]
                coeff = coeff[:,min_idx]
                
                coeff_strength = np.sum(np.abs(coeff))
                coeff_number = np.sum(coeff != 0)
                best_score = - coeff_strength - coeff_number
                
                coef[run_number,:] = np.transpose(coeff)
                seed_vec[run_number] = cache.seed
                stop_condition = 0
                
            elif conflict_new[min_idx] == best_conflict:
                
                coeffs = new_s[:,condition]
                coeffs = coeffs[:, conflict_new == conflict_new[min_idx] ]
                
                coeffs_strength = np.sum(np.abs(coeffs), 0)
                coeffs_number = np.sum(coeffs != 0, 0)
                score = - coeffs_strength - coeffs_number
                
                max_idx = np.argmax(score)
                
                if score[max_idx] > best_score:
                    best_score = score[max_idx]
                    coef[run_number,:] = np.transpose( coeffs[:,max_idx] )
                    seed_vec[run_number] = cache.seed
                    stop_condition = 0
        
        states.append(np.transpose(s))
        actions.append(a)
        total_reward += r
        
        s = new_s
    
    stop_condition += 1
    flag = True if stop_condition == 40 else False
    
    return states, actions, total_reward, flag, best_score, best_conflict

def select_elites(states_batch, actions_batch, rewards_batch, coeffs_size, percentile=50):
    """
    Select states and actions from games that have rewards >= percentile
    :param states_batch: list of lists of states, states_batch[session_i][t]
    :param actions_batch: list of lists of actions, actions_batch[session_i][t]
    :param rewards_batch: list of rewards, rewards_batch[session_i][t]
    
    :returns: elite_states,elite_actions, both 1D lists of states and respective actions from elite sessions
    """
    
    reward_threshold = np.percentile(rewards_batch, percentile)
    mask = (rewards_batch >= reward_threshold)
    
    states_batch = np.array(states_batch).transpose(1,0,2).reshape(n_sessions,-1)
    elite_states = states_batch[mask].reshape(-1,coeffs_size)
    
    actions_batch = np.transpose(actions_batch)
    elite_actions = actions_batch[mask].reshape(-1,1)
    
    return elite_states, elite_actions

def show_progress(batch_rewards, percentile):
    
    mean_reward, threshold = np.mean(batch_rewards), np.percentile(batch_rewards, percentile)
    #print("mean reward = %.3f, threshold=%.3f"%(mean_reward, threshold))
    #Log.append([int(mean_reward), int(threshold)])
    
    print("mean reward = %.1f, threshold=%.1f"%(mean_reward, threshold))
    return mean_reward
    '''
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
    '''


############################# MAIN FUNCTION ###################################


def get_envelopes(pretrain_idx = -1, given_seed = []):
    global  stop_condition, coef, found, run_number, seed_vec
    start_runtime = time.time()
    
    with open(filein,"r") as f:
        data = f.readlines()
    num_cases = int(len(data)/2)
    for Case in range(num_cases):
        if Case:
            pretrain_idx = -1
            given_seed = []
        start_time = time.time()
    
        LHS_string = data[2*Case]
        aux = int( data[2*Case + 1] )
        
        [init_training, LHS, allbits, reset_state, n, pretrain_idx] = pretrain(aux, LHS_string, pretrain_idx)
        
        coeffs_size = int( n*(n+1)/2 + 1)
        n_actions = 2*coeffs_size + 1
        
        LHS_size = 2**(n-aux)
        
        LHS.shape = (LHS_size, 1)
        reset_state.shape = (coeffs_size, 1)
        reset_state = np.repeat(reset_state, n_sessions, axis = 1)
        
        one_hot = np.zeros((coeffs_size,n_actions))
        for i in range(coeffs_size):
            one_hot[i,i] = 1
            one_hot[i,coeffs_size+i] = -1
        
        env = environment(one_hot)
        
        cache = create_cache(LHS, allbits, coeffs_size, n_actions, n, aux, reset_state, env)
    
        # create agent
        
        agent = MLPClassifier(
            hidden_layer_sizes = (coeffs_size),
            activation='tanh',
            warm_start=False,#True,  # keep progress between .fit() calls
            max_iter=1,       # make only 1 iteration on each .fit()
        )
        
        init_training = np.array(init_training)
        size = init_training.shape[0]
        
        train_data = np.subtract(init_training,one_hot[:,0])
        target = np.ones((size,1))*0
        for action in range(1,n_actions):
            train_data = np.append( train_data, np.subtract(init_training,one_hot[:,action]) )
            target = np.append( target, np.ones((size,1))*action )
        
        train_data = train_data.reshape(-1,coeffs_size)
        
        # initialize agent
        #pretrained_agent = agent.fit(train_data,target)
        
        #######################################################################
        # Finding COMIGs
        
        best_num_of_comigs = 1000
        for Run in range(1):
            coef = np.zeros((max_runs,coeffs_size))
            seed_vec = np.zeros(max_runs)
            score = 0
            run_number = 0
            while True:
                
                mask = np.array([True]*LHS_size)
                mask = mask.reshape(-1,1)
                mask_preserve = np.zeros((max_runs,LHS_size,1))
                
                for i in np.arange(run_number):
                    RHS = rhs(coef[i].reshape(-1,), cache.allbits, cache.aux).reshape(-1,1)
                    mask_preserve[i] = np.subtract(RHS,LHS) != 0
                    mask = np.multiply(mask, mask_preserve[i])
                
                mask = mask.reshape(-1,)
                mask = mask == 1
                mask_sum = np.sum(mask)
                if  not mask_sum or run_number == max_runs:
                    break
                print('states left = ', mask_sum)
                
                cache.add_mask(mask, mask_sum)
                
                best_run_score = -10000
                best_conflict = 1
                for _ in range(4):
                    if given_seed:
                        seed = given_seed[run_number]
                    else:
                        seed = np.random.randint(1000000)
                    np.random.seed(seed)
                    cache.add_seed(seed)
                    agent.fit(train_data,target)
                    for j in range(100):
                        #mean_reward_max_in = -10000
                        
                        sessions = generate_session(agent, cache, best_run_score, best_conflict, t_max)
                        
                        batch_states, batch_actions, batch_rewards, flag, best_run_score, best_conflict = sessions[0], sessions[1], sessions[2], sessions[3], sessions[4], sessions[5]
                        
                        if flag:
                            break
                        
                        elite_states, elite_actions = select_elites(batch_states, batch_actions, batch_rewards, cache.coeffs_size, percentile)
                        
                        #mean_reward = show_progress(batch_rewards, percentile)
                        
                        agent.partial_fit(elite_states, elite_actions)
                        
                    if given_seed:
                        break
                
                if best_run_score > -10000:
                    run_number += 1
                    print('Case: ',Case,'Run: ',Run,'Comig: ',run_number)
                    score = score + best_run_score
            print('num of comigs = ',run_number)
            if run_number < best_num_of_comigs:
                best_num_of_comigs = run_number
                Best_score = score
                best_coef = coef[:run_number]
                #savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
            elif run_number == best_num_of_comigs:
                if score > Best_score:
                    Best_score = score
                    best_coef = coef[:run_number]
                    #savemat('data.mat',{ ('coef' + str(Case)) :best_coef} )
        
        print('Done!')
        print('best num of comigs = ',best_num_of_comigs,'\n')
        
        # verify coef and get const_terms
        #[verified, const_terms] = eng.verify(n, aux, LHS.tolist(), allbits.tolist(), best_coef.tolist(), nargout = 2)
        
        with open(fileout + str(Case) + ".txt","a+") as file1:
            
            if (np.min(rhs(np.transpose(best_coef), cache.allbits, cache.aux),1) != np.transpose(LHS)).any():
                file1.write('\nNot verified!!!\n')
            
            file1.write(LHS_string[:-1] + ':\n\n')
            
            max_percent = -1
            max_idx = -1
            for m in range(best_coef.shape[0]):
                mask = ( rhs(best_coef[m], cache.allbits, cache.aux) != np.transpose(LHS) )
                percent = 100 * ( 1 - np.sum(np.sum(mask))/LHS.shape[0] )
                if percent > max_percent:
                    max_percent = percent
                    max_idx = m
            
            if max_idx:
                print('First run is NOT the best')
                print('Best run: %d, with percentage: %.2f', max_idx, max_percent)
            
            mask = np.ones((1,LHS.shape[0]))            
            for m in range(best_coef.shape[0]):
                mask = mask * ( rhs(best_coef[m], cache.allbits, cache.aux) != np.transpose(LHS) )
                percent = 100 * ( 1 - np.sum(np.sum(mask))/LHS.shape[0] )
                
                file1.write('(%6.2f%%)  ' % percent )
                vec = best_coef[m,:]
                k = 0
                for i in range(1,n):
                    for j in range(i+1,n+1):
                        if vec[k] != 0:
                            if vec[k] == 1:
                                file1.write(' + b_{%d}b_{%d}' % (i, j) )
                            elif vec[k] == -1:
                                file1.write(' - b_{%d}b_{%d}' % (i, j) )
                            else:
                                file1.write(' %+d b_{%d}b_{%d}' % (vec[k], i, j) )
                        k += 1
                
                for i in range(1,n+1):
                    if vec[k] != 0:
                        if vec[k] == 1:
                            file1.write(' + b_{%d}' % i)
                        elif vec[k] == -1:
                            file1.write(' - b_{%d}' % i)
                        else:
                            file1.write(' %+d b_{%d}' % (vec[k], i) )
                    k += 1
                
                if vec[k] != 0:
                    file1.write(' %+d' % vec[k] )
                file1.write('\n')
            
            file1.write('\nTotal Runtime: %6.2f seconds' % (time.time() - start_time) )
            
            file1.write('\n\nSeeds used:\nPretrain index: %6d\n' % pretrain_idx )
            for idx in range(best_num_of_comigs):
                file1.write('         Run %d: %6d\n' % (idx + 1, seed_vec[idx]) )
            file1.write('\n')
    
    total = time.time() - start_runtime
    print('total runtime: %.1f minutes' % float(total/60) )
    print('avg runtime per comig = %.1f minutes' % float(total/60/num_cases) )


if profiling == 0:
    get_envelopes()

elif profiling == 1:
    import cProfile, pstats, io
    
    pr = cProfile.Profile()
    pr.enable()
    get_envelopes()
    pr.disable()
    
    s = io.StringIO()
    sortby = 'cumulative'
    
    temp = pstats.Stats(pr)
    temp.dump_stats("output.pstats")
    
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats(20)
    print(s.getvalue())

elif profiling == 2:
    import line_profiler
    
    l = line_profiler.LineProfiler()
    l.add_function(get_envelopes)
    l.add_function(generate_session)
    l.add_function(rhs)
    l.add_function(reward)
    l.add_function(accuracy)
    l.add_function(dec2base)
    l.add_function(pretrain)
    l.run('get_envelopes()')
    l.print_stats()


