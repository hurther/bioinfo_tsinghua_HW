import numpy as np
from hmmlearn import hmm
import sys
import random

sys.path.append("./src")

from HMMv1 import HMMv1

# fix random seed
np.random.seed(42)
random.seed(42)

states = ["box 1", "box 2", "box3"]
n_states = len(states)

observations = ["red", "white"]
n_observations = len(observations)

ip = np.array([0.2, 0.4, 0.4])
A = np.array([[0.5, 0.2, 0.3], 
                [0.3, 0.5, 0.2],
                [0.2, 0.3, 0.5]])
B = np.array([[0.5, 0.5],
                [0.4, 0.6],
                [0.7, 0.3]])

def init_hmm(nStates, nObs):
    """
    Initialize HMM model.
    
    `nStates`: number of states
    
    `nObs`: number of observe
    
    return: ip, a, b
    """
    ip = np.random.rand(nStates)
    ip = ip / np.sum(ip)
    a = np.random.rand(nStates, nStates)
    a = a / np.sum(a, axis=1, keepdims=True)
    b = np.random.rand(nStates, nObs)
    b = b / np.sum(b, axis=1, keepdims=True)
    return ip, a, b

# random init
ip_random, a_random, b_random = init_hmm(n_states, n_observations)
print("ip_random:", ip_random)
print("a_random:", a_random)
print("b_random:", b_random)

print("-"*20, "hmmlearn", "-"*20)

model = hmm.MultinomialHMM(n_components=n_states, n_iter=20, tol=0.01)
model.startprob_ = ip
model.transmat_ = A
model.emissionprob_ = B

# predict a sequence of hidden states based on visible states
seen = np.array([[0,1,0]])
logprob, box = model.decode(seen, algorithm="viterbi")
print("The ball picked:", ", ".join(map(lambda x: observations[x], seen.flatten())))
print("The hidden box:", ", ".join(map(lambda x: states[x], box)))

# predict the probability of a sequence of visible states
print(model.score(seen))
print(np.exp(model.score(seen)))

# fit model parameters based on data
X = np.array([[0,1,0,1], [0,0,0,1], [1,0,1,1]])

new_model = hmm.MultinomialHMM(n_components=n_states, n_iter=100, tol=0.01, init_params=" ")
new_model.startprob_ = ip_random
new_model.transmat_ = a_random
new_model.emissionprob_ = b_random
new_model.n_features = n_observations
new_model.fit(X, lengths=[1, 1, 1]) # if not set lengths, X will be considered as a single sequence
print(new_model.startprob_)
print(new_model.transmat_)
print(new_model.emissionprob_)
print(new_model.score(X[0].reshape(1, -1)))

print("-"*20, "HMMv1", "-"*20)

model = HMMv1(init_probability=ip, transition_probability=A, emission_probability=B, n_iter=20, tol=0.01)

# predict a sequence of hidden states based on visible states
seen = np.array([0,1,0])
box = model.decode(seen)
print("The ball picked:", ", ".join(map(lambda x: observations[x], seen)))
print("The hidden box", ", ".join(map(lambda x: states[x], box)))

# predict the probability of a sequence of visible states
print(np.log(model.score(seen)))
print(model.score(seen))

# fit model parameters based on data
new_model = HMMv1(init_probability=ip_random, transition_probability=a_random, emission_probability=b_random, n_iter=100, tol=0.01)
new_model.fit(X)
print(new_model.ip)
print(new_model.a)
print(new_model.b)
print(np.log(new_model.score(X[0])))