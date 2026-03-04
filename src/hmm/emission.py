"""
Emission model for Hidden Markov Model ancestry inference.

Calculates the probability of observing a genotype given the ancestry state (YRI or CEU).
Uses Hardy-Weinberg equilibrium: P(g|p) where p is allele frequency and g is genotype.
"""

import numpy as np

class EmissionModel:
    def __init__(self, yri_freqs, ceu_freqs, epsilon=1e-6):
        """
        yri_freqs: dict {pos: freq}
        ceu_freqs: dict {pos: freq}
        epsilon: small floor value to avoid zero probabilities
        """
        self.yri_freqs = yri_freqs
        self.ceu_freqs = ceu_freqs
        self.epsilon = epsilon

    def get_log_emission(self, state, genotype, position):
        # Get the frequency for the given state (YRI or CEU)
        freqs = self.yri_freqs if state == "YRI" else self.ceu_freqs
        
        # If the SNP position isn't in our frequency dict, return a neutral log-prob
        if position not in freqs:
            return np.log(0.5)
            
        p = freqs[position]
        
        # Basic emission logic: 
        # Probability of observing genotype (0,0), (0,1), or (1,1) given allele frequency p
        # We add epsilon to prevent log(0) errors
        if genotype == (0, 0):
            prob = (1 - p)**2
        elif genotype == (0, 1) or genotype == (1, 0):
            prob = 2 * p * (1 - p)
        elif genotype == (1, 1):
            prob = p**2
        else:
            prob = 0.5 # Fallback
            
        return np.log(max(prob, self.epsilon))

    def get_emission_probs(self, pos, genotype):
        """
        Calculates P(Genotype | State) for both YRI and CEU.
        genotype: tuple (0, 1), (1, 1), etc.
        """
        f_yri = self.yri_freqs.get(pos, 0.5) # Default to 0.5 if SNP missing
        f_ceu = self.ceu_freqs.get(pos, 0.5)
        
        # Calculate likelihood for YRI state
        # We use a small epsilon so we never multiply by exactly 0
        p_yri = self._calc_likelihood(genotype, f_yri)
        
        # Calculate likelihood for CEU state
        p_ceu = self._calc_likelihood(genotype, f_ceu)
        
        return {"YRI": p_yri, "CEU": p_ceu}

    def _calc_likelihood(self, genotype, f):
        """
        Diploid likelihood: P(a1|f) * P(a2|f)
        """
        prob = 1.0
        for allele in genotype:
            if allele == 1:
                prob *= f
            else:
                prob *= (1 - f)
        return max(prob, self.epsilon)