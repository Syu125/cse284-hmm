"""
Transition model for Hidden Markov Model ancestry inference.

Calculates transition probabilities between ancestry states based on genetic distance.
Uses recombination rate to estimate probability of switching states between SNPs.
"""

import numpy as np

class TransitionModel:
    def __init__(self, generations=10):
        self.G = generations # Generations since admixture
        # Haplotype-level ancestry states (K states, K=2 for this project).
        self.states = ["CEU", "YRI"]

    def get_transition_matrix(self, cm_start, cm_end):
        """
        Calculates transition probabilities among haplotype ancestry states.
        cm_start: genetic position of SNP i
        cm_end: genetic position of SNP i+1
        """
        # 1. Calculate distance in Morgans (cM / 100)
        d = abs(cm_end - cm_start) / 100.0
        
        # 2. Probability of a switch
        # If distance is 0, p_switch is 0.
        p_switch = 1 - np.exp(-self.G * d)
        
        p_stay = 1 - p_switch

        # State order: CEU, YRI
        matrix = np.array(
            [
                [p_stay, p_switch],
                [p_switch, p_stay],
            ],
            dtype=float,
        )

        return matrix