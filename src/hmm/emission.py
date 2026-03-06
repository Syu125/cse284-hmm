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
        Calculates P(Genotype | State) for diploid ancestry states.
        genotype: tuple (0, 1), (1, 1), etc.
        """
        f_yri = self.yri_freqs.get(pos, 0.5) # Default to 0.5 if SNP missing
        f_ceu = self.ceu_freqs.get(pos, 0.5)
        
        p_ceu_ceu = self._calc_likelihood(genotype, f_ceu, f_ceu)
        p_yri_yri = self._calc_likelihood(genotype, f_yri, f_yri)
        p_ceu_yri = self._calc_mixed_likelihood(genotype, f_ceu, f_yri)

        return {
            "CEU_CEU": p_ceu_ceu,
            "CEU_YRI": p_ceu_yri,
            "YRI_YRI": p_yri_yri,
        }

    def _calc_likelihood(self, genotype, f1, f2):
        """
        Diploid likelihood with potentially different ancestral frequencies per allele.
        """
        if genotype is None or len(genotype) < 2:
            return 0.5

        allele_1 = genotype[0]
        allele_2 = genotype[1]

        p1 = self._allele_prob(allele_1, f1)
        p2 = self._allele_prob(allele_2, f2)
        prob = p1 * p2
        return max(prob, self.epsilon)

    def _calc_mixed_likelihood(self, genotype, f_left, f_right):
        """Average over the two possible allele-origin orderings for mixed ancestry."""
        prob_a = self._calc_likelihood(genotype, f_left, f_right)
        prob_b = self._calc_likelihood(genotype, f_right, f_left)
        return max(0.5 * (prob_a + prob_b), self.epsilon)

    def _allele_prob(self, allele, f):
        if allele == 1:
            return f
        if allele == 0:
            return 1 - f
        return 0.5