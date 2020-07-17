
# ###### An estimate of basic secondary structure in protein (amino acid) sequences. The model we consider is a simplistic rendition of the model discussed in S C. Schmidler et al. (2004) Bayesian Segmentation of Protein Secondary Structure, doi:10.1089/10665270050081496
#https://github.com/adriansunachoi/phylogeneticTreeReconstruction/



get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import random, math
from Node import Node
from Tree import Tree
import tree69
states = np.array(['H','S','T'])
symbols = np.array(['N','B','I'])
emission = np.array([[0.1,0.3,0.6],[0.3,0.55,0.15],[0.7,0.1,0.2]])
transition = np.array([[15/16,3/160,7/160],[2/45,8/9,1/15],[1/14,1 /14,6/7]])





def sequence_simulate(n=150):
    '''
    To simulate state & symbol sequences of arbirary length from the HMM
    To print out a simulated sequence each for state & symbol

    :param <n>: sequence length
    :return <state_sequence>, <symbol_sequence>, &type=numpy arrays

    '''
    symbol_sequence = []
    
    # The first state is chosen at an equal chance of one third
    state_sequence = [states[random.choice([0,1,2])]]
    for i in range(n):
        
        # A state in HMM is determined by only the previous state
        # The last element in the state sequence is used
        current_state = state_sequence[-1]
        
        # Find the index of the current index {0,1,2}
        current_state_index = int(np.where(states==current_state)[0])
        
        # A symbol is emitted based on the respective emission probability
        current_symbol = np.random.choice(symbols, 1, p=emission[current_state_index])[0]
        
        # The emitted symbol is appended to the sequence
        symbol_sequence.append(current_symbol)
        
        # Transition to next state based on the respective transition probability
        next_state = np.random.choice(states, 1, p=transition[current_state_index])[0]
        
        # The new state is appended to the sequence
        # Since a new state is always added to the sequence, when reach the required length, break the loop
        if i==n-1:
            break
        state_sequence.append(next_state)
        
    return state_sequence, symbol_sequence





#####TEST#####
test_state_seq, test_symbol_seq = sequence_simulate()
print('State sequence for length=150\n',''.join(test_state_seq))
print('Symbol sequence for length=150\n',''.join(test_symbol_seq))


# ##### Natural log of joint probability is calculated



def joint_prob(pi,x):
    '''
    To calcultae the natural log of the joint probability P(pi, x)
    
    : param <pi>: a state sequence
            <x>: a symbol sequence
            &dtype=list
    : return <math.log(result)>: the natural log of the joint probability
    
    '''
    # Initial prob of choosing any state is equally one third
    initial_prob = result = 1/3 
    
    # Enumerate through all states and symbols
    for i in range(0, len(x)-1):
        current_state_index = int(np.where(states==pi[i])[0])
        next_state_index = int(np.where(states==pi[i+1])[0])
        current_symbol_index = int(np.where(symbols==x[i])[0])
        
        # Multiply all terms
        result*=emission[current_state_index, current_symbol_index]*transition[current_state_index, next_state_index]
    return math.log(result)




S = 'S'
H = 'H'
T = 'T'
B = 'B'
I = 'I'
N = 'N'
pi = [S,S,H,H,H,T,T,S,S,S,H,H,H,H,H,H,S,S,S,S,S,S]
x = [B,I,N,B,N,I,N,B,N,I,N,B,I,N,B,I,I,N,B,B,N,B]

print("Joint probability for pi & x\n",joint_prob(pi,x))
print("Joint probability for simulated sequences in b\n", joint_prob(test_state_seq, test_symbol_seq))

