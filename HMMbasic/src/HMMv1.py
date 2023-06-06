import numpy as np

class HMMv1:
    """    
    A simple implementation of HMM with discrete observations.
    
    Attention:
        
        If the observations is too long, it may cause underflow.
        
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
        
        Attention: if the observations is too long, it may cause underflow.
        
        Args:
            obs: A list of observations.
            
        Returns:
            alpha: A matrix of forward probabilities.
            
            prob: The probability of the observations.
        """
        nStates = np.shape(self.b)[0]
        T = np.shape(obs)[0]
        alpha = np.zeros((nStates, T))
        alpha[:, 0] = self.ip * self.b[:, obs[0]]
    
        for t in range(1, T):
            for s in range(nStates):
                alpha[s, t] = self.b[s, obs[t]] * np.sum(alpha[:, t-1] * self.a[:, s])
            
        return alpha, np.sum(alpha[:, -1])
        
    def backward(self, obs=None):
        """
        Backward algorithm.
        
        Attention: if the observations is too long, it may cause underflow.
        
        Args:
            obs: A list of observations.
            
        Returns:
            beta: A matrix of backward probabilities.
        
            prob: The probability of the observations.
        """
        nStates = np.shape(self.b)[0]
        T = np.shape(obs)[0]
        beta = np.zeros((nStates, T))
        beta[:, -1] = 1.0
        
        for t in range(T-2, -1, -1):
            for s in range(nStates):
                beta[s, t] = np.sum(beta[:, t+1] * self.a[s, :] * self.b[:, obs[t+1]])
                
        return beta, np.sum(self.ip * self.b[:, obs[0]] * beta[:, 0])
    
    def score(self, obs=None):
        """
        Calculate the probability of the observations.
        
        Args:
            obs: A list of observations.
            
        Returns:
            prob: The probability of the observations.
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
        delta = np.zeros((nStates, T))
        psi = np.zeros((nStates, T))
        delta[:, 0] = self.ip * self.b[:, obs[0]]
        psi[:, 0] = 0
        
        for t in range(1, T):
            for s in range(nStates):
                delta[s, t] = np.max(delta[:, t-1] * self.a[:, s]) * self.b[s, obs[t]]
                psi[s, t] = np.argmax(delta[:, t-1] * self.a[:, s])
                
        path = np.zeros(T, dtype=np.int64)
        path[-1] = np.argmax(delta[:, -1])
        
        for t in range(T-2, -1, -1):
            path[t] = psi[int(path[t+1]), t+1]
            
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
        nStates = np.shape(self.b)[0]
        nObservations = np.shape(self.b)[1]
        nSeq = len(X)
        nIts = 0
        while True:
            nIts += 1
            old_a = self.a.copy()
            old_b = self.b.copy()
            T = np.shape(X)[1]
            total_xi = np.zeros((nStates, nStates, T))
            total_gamma = np.zeros((nStates, T))
            
            for seq in X:
                xi = np.zeros((nStates, nStates, T))
                alpha = self.forward(seq)[0]
                beta = self.backward(seq)[0]
                gamma = alpha * beta
                gamma /= gamma.sum(0)
                
                # E-step
                for t in range(T-1):
                    for i in range(nStates):
                        for j in range(nStates):
                            xi[i, j, t] = alpha[i, t] * self.a[i, j] * self.b[j, seq[t+1]] * beta[j, t+1]
                    xi[:, :, t] /= xi[:, :, t].sum()
                    
                # the last step has no beta, b
                for i in range(nStates):
                    for j in range(nStates):
                        xi[i, j, -1] = alpha[i, -1] * self.a[i, j]
                xi[:, :, -1] /= xi[:, :, -1].sum()
                
                total_xi += xi
                total_gamma += gamma
            
            # M-step
            for i in range(nStates):
                for j in range(nStates):
                    self.a[i, j] = np.sum(total_xi[i, j, :]) / np.sum(total_gamma[i, :-1])
                    
            self.a /= self.a.sum(1, keepdims=True)
            
            for j in range(nStates):
                for k in range(nObservations):
                    found = (seq == k).nonzero()[0]
                    self.b[j, k] = np.sum(total_gamma[j, found]) / np.sum(total_gamma[j])
                    
            self.b /= self.b.sum(1, keepdims=True)
            self.ip = total_gamma[:, 0] / nSeq
            
            # check convergence
            if np.linalg.norm(self.a - old_a) < self.tol and np.linalg.norm(self.b - old_b) < self.tol:
                print('Converged after %d iterations.' % nIts)
                break
            
            # check max iterations
            if nIts >= self.n_iter:
                print('Max iterations reached.')
                break
        return self