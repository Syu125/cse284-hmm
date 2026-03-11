"""
Transition model for Hidden Markov Model ancestry inference.

Calculates transition probabilities between ancestry states based on genetic distance.
Uses recombination rate to estimate probability of switching states between SNPs.
"""

import numpy as np

class TransitionModel:
    def __init__(self, generations=10):
        self.G = generations # Generations since admixture

    def get_transition_matrix(self, cm_start, cm_end):
        """
        Calculates transition probabilities among phased diploid ancestry states.
        cm_start: genetic position of SNP i
        cm_end: genetic position of SNP i+1
        """
        # 1. Calculate distance in Morgans (cM / 100)
        d = abs(cm_end - cm_start) / 100.0
        
        # 2. Probability of a switch
        # If distance is 0, p_switch is 0.
        p_switch = 1 - np.exp(-self.G * d)
        
        p_stay = 1 - p_switch

        # State order: CEU_CEU, CEU_YRI, YRI_CEU, YRI_YRI
        haplotypes = ["CEU", "YRI"]
        states = [
            ("CEU", "CEU"),
            ("CEU", "YRI"),
            ("YRI", "CEU"),
            ("YRI", "YRI"),
        ]

        hap_transition = {
            (src, dst): (p_stay if src == dst else p_switch)
            for src in haplotypes
            for dst in haplotypes
        }

        matrix = np.zeros((len(states), len(states)), dtype=float)
        for i, (h1_src, h2_src) in enumerate(states):
            for j, (h1_dst, h2_dst) in enumerate(states):
                matrix[i, j] = hap_transition[(h1_src, h1_dst)] * hap_transition[(h2_src, h2_dst)]

        return matrix