import numpy as np
from hmmlearn import hmm
import sys
import random

sys.path.append("./src")

from HMMv1 import HMMv1
from HMMv2 import HMMv2

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

O = np.array([0,1,0]).repeat(1000) # the observation sequence is 3000 long

# hmmlearn
print('-'*20 + 'hmmlearn' + '-'*20)
model = hmm.MultinomialHMM(n_components=n_states)
model.startprob_ = ip
model.transmat_ = A
model.emissionprob_ = B

print(model.score(O.reshape(1,-1)))

# HMMv1
print('-'*20 + 'HMMv1' + '-'*20)
model = HMMv1(ip, A, B)

print(model.score(O))

# HMMv2
print('-'*20 + 'HMMv2' + '-'*20)
model = HMMv2(ip, A, B)

print(model.score(O))