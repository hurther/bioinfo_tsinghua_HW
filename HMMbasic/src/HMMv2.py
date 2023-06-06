import numpy as np
from utils import *

class HMMv2:
    """
    A simple implementation of HMM with discrete observations.
    
    Considering the overflow problem, we use log probability instead of probability.
    
    """
    def __init__(self, init_probability=None, transition_probability=None, emission_probability=None, n_iter=100, tol=1e-2):
        """
        Args:
            observe_sequences: A list of observations.
            
            init_probability: A list of initial probabilities.
            
            transition_probability: A matrix of transition probabilities.
            
            emission_probability: A matrix of emission probabilities.
            
            n_iter: The maximum number of iterations.
            
            tol: The tolerance of the difference between the log likelihood of the current model and the previous model.
        """
        self.n_iter = n_iter
        self.tol = tol
        self.ip = init_probability
        self.a = transition_probability
        self.b = emission_probability
        
    def forward(self, obs=None):
        """
        Forward algorithm.
        
        Args:
            obs: A list of observations.
            
        Returns:
            alpha: A matrix of forward probabilities.
            
            prob: The probability of the observations in log.
        """
        nStates = np.shape(self.b)[0]
        T = np.shape(obs)[0]
        
        work_buffer = np.zeros(nStates)
        framelogprob = log_mask_zero(self.b)[:,obs].T
        log_init_prob = log_mask_zero(self.ip)
        log_transmat = log_mask_zero(self.a)
        alpha = np.zeros((T, nStates))
        
        for i in range(nStates):
            alpha[0, i] = log_init_prob[i] + framelogprob[0, i]
            
        for t in range(1, T):
            for j in range(nStates):
                for i in range(nStates):
                    work_buffer[i] = alpha[t-1, i] + log_transmat[i, j]
                alpha[t, j] = logsumexp(work_buffer) + framelogprob[t, j]
                
        return alpha, logsumexp(alpha[-1])
    
    def backward(self, obs=None):
        """
        Backward algorithm.
        
        Args:
            obs: A list of observations.
            
        Returns:
            beta: A matrix of backward probabilities.
        
            prob: The probability of the observations in log.
        """
        nStates = np.shape(self.b)[0]
        T = np.shape(obs)[0]
        
        work_buffer = np.zeros(nStates)
        framelogprob = log_mask_zero(self.b)[:,obs].T
        log_init_prob = log_mask_zero(self.ip)
        log_transmat = log_mask_zero(self.a)
        beta = np.zeros((T, nStates))
        
        for i in range(nStates):
            beta[nStates-1, i] = 0.0
            
        for t in range(T-2, -1, -1):
            for i in range(nStates):
                for j in range(nStates):
                    work_buffer[j] = log_transmat[i, j] + framelogprob[t+1, j] + beta[t+1, j]
                beta[t, i] = logsumexp(work_buffer)
                
        return beta, logsumexp(log_init_prob + framelogprob[0, :] + beta[0, :])

    def score(self, obs=None):
        """
        Calculate the probability of the observations.
        
        Args:
            obs: A list of observations.
            
        Returns:
            prob: The probability of the observations in log.
        """
        _, prob = self.forward(obs)
        return prob

    def decode(self, obs=None):
        """
        Viterbi algorithm.
        
        Args:
            obs: A list of observations.
            
        Returns:
            path: The most possible path.
        """
        nStates = np.shape(self.b)[0]
        T = np.shape(obs)[0]
        
        path = np.empty(T, dtype=int)
        work_buffer = np.zeros(nStates)
        framelogprob = log_mask_zero(self.b)[:,obs].T
        log_init_prob = log_mask_zero(self.ip)
        log_transmat = log_mask_zero(self.a)
        delta = np.empty((T, nStates))
        
        for i in range(nStates):
            delta[0, i] = log_init_prob[i] + framelogprob[0, i]
            
        for t in range(1, T):
            for j in range(nStates):
                for i in range(nStates):
                    work_buffer[i] = delta[t-1, i] + log_transmat[i, j]
                delta[t, j] = np.max(work_buffer) + framelogprob[t, j]
                
        psi = path[T-1] = np.argmax(delta[T-1])
        
        for t in range(T-2, -1, -1):
            for i in range(nStates):
                work_buffer[i] = delta[t, i] + log_transmat[i, psi]
            psi = path[t] = np.argmax(work_buffer)
            
        return path
    
    def init_hmm(self, nStates=None, nObservations=None):
        """
        Initialize the parameters of HMM.
        
        Args:
            nStates: The number of hidden states.
            
            nObservations: The number of observations.
        """
        if nStates is None:
            nStates = np.shape(self.a)[0]
        if nObservations is None:
            nObservations = np.shape(self.b)[1]
        ip = np.random.rand(nStates)
        self.ip = ip / ip.sum()
        a = np.random.rand(nStates, nStates)
        self.a = a /np.sum(a, axis=1, keepdims=True)
        b = np.random.rand(nStates, nObservations)
        self.b = b / np.sum(b, axis=1, keepdims=True)
        return self

    
    def fit(self, X):
        """
        Baum-Welch algorithm for fitting HMM with the observations.
        
        Args:
            X: A 2-D list of observations. [nSeq, T]
        """
        nFeatures = np.shape(self.b)[1]
        nStates = np.shape(self.b)[0]
        history = []
        for iter in range(self.n_iter):
            startprob = np.zeros(nStates)
            trans = np.zeros((nStates, nStates))
            emission = np.zeros((nStates, nFeatures))
            curr_logprob = 0
            # E-step
            for seq in X:
                framelogprob = log_mask_zero(self.b)[:,seq].T
                alpha, logprob = self.forward(seq)
                curr_logprob += logprob
                beta, _ = self.backward(seq)
                
                gamma = alpha + beta
                log_normalize(gamma, axis=1)
                with np.errstate(under='ignore'):
                    posteriors = np.exp(gamma)
                    
                n_samples, n_components = framelogprob.shape
                startprob += posteriors[0]
                xi_sum = np.full((n_components, n_components), -np.inf)
                work_buffer = np.full((n_components, n_components), -np.inf)
                logprob = logsumexp(alpha[n_samples-1])
                log_transmat = log_mask_zero(self.a)
                
                for t in range(n_samples-1):
                    for i in range(n_components):
                        for j in range(n_components):
                            work_buffer[i, j] = alpha[t, i] + log_transmat[i, j] + framelogprob[t+1, j] + beta[t+1, j] - logprob
                            
                    for i in range(n_components):
                        for j in range(n_components):
                            xi_sum[i, j] = logsumexp([xi_sum[i, j], work_buffer[i, j]])
                with np.errstate(under='ignore'):
                    trans += np.exp(xi_sum)
                    
                for t, symbol in enumerate(seq):
                    emission[:, symbol] += posteriors[t]
            # M-step
            startprob_ = np.maximum(startprob, 0)
            self.ip = np.where(self.ip == 0, 0, startprob_)
            normalize(self.ip)
            
            transmat_ = np.maximum(trans, 0)
            self.a = np.where(self.a == 0, 0, transmat_)
            normalize(self.a, axis=1)
            
            self.b = emission / np.sum(emission, axis=1, keepdims=True)
        
            history.append(curr_logprob)
            # check for convergence
            if (len(history) >= 2 and history[-1] - history[-2] < self.tol and history[-2] - history[-3] < self.tol):
                print('Converged at iteration %d.' % iter)
                break
            if iter == self.n_iter-1:
                print('Max iteration reached.')
                break
        
        return self