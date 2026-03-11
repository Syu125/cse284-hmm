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
        # Haplotype-level ancestry states (K states, K=2 for this project).
        self.states = ["CEU", "YRI"]

    def get_haplotype_emission_probs(self, pos, allele):
        """
        Calculates P(allele | ancestry_state) for a single phased haplotype.
        allele: 0 or 1.
        """
        f_yri = self.yri_freqs.get(pos, 0.5)
        f_ceu = self.ceu_freqs.get(pos, 0.5)

        return {
            "CEU": max(self._allele_prob(allele, f_ceu), self.epsilon),
            "YRI": max(self._allele_prob(allele, f_yri), self.epsilon),
        }

    def get_emission_probs(self, pos, genotype):
        """
        Backward-compatible diploid ancestry likelihoods for homogeneous states.
        genotype: tuple (0, 1), (1, 1), etc.

        Returns:
            dict with keys {"CEU", "YRI"}.
        """
        f_yri = self.yri_freqs.get(pos, 0.5) # Default to 0.5 if SNP missing
        f_ceu = self.ceu_freqs.get(pos, 0.5)

        return {
            "CEU": self._calc_likelihood(genotype, f_ceu, f_ceu),
            "YRI": self._calc_likelihood(genotype, f_yri, f_yri),
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

    def _allele_prob(self, allele, f):
        if allele == 1:
            return f
        if allele == 0:
            return 1 - f
        return 0.5