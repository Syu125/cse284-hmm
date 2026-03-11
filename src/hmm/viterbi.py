"""
Viterbi algorithm for Hidden Markov Model ancestry inference.

Implements the Viterbi algorithm to find the most likely sequence of ancestry states
for a given sequence of genotypes. Uses dynamic programming to efficiently compute
the maximum probability path through the HMM state space.
"""

import numpy as np

class InferenceEngine:
    def __init__(self, emission_model, transition_model):
        self.emissions = emission_model
        self.transitions = transition_model
        self.states = ["CEU_CEU", "CEU_YRI", "YRI_CEU", "YRI_YRI"]

    def run_viterbi(self, snp_positions, genotypes, genetic_map_func):
        """
        snp_positions: list of physical positions
        genotypes: list of tuples (0,1), (1,1), etc.
        genetic_map_func: a lambda or function that takes pos and returns cM
        """
        n_snps = len(snp_positions)
        if n_snps == 0:
            return []

        n_states = len(self.states)
        
        # viterbi[state, snp] stores the max probability to reach that point
        viterbi = np.zeros((n_states, n_snps))
        # backpointer stores which state we came from
        backpointer = np.zeros((n_states, n_snps), dtype=int)

        # 1. Initialize (First SNP)
        first_pos = snp_positions[0]
        first_em = self.emissions.get_emission_probs(first_pos, genotypes[0])
        initial_prior = 1.0 / n_states
        for state_idx, state_name in enumerate(self.states):
            viterbi[state_idx, 0] = np.log(initial_prior * first_em[state_name])

        # 2. Recursion
        for i in range(1, n_snps):
            curr_pos = snp_positions[i]
            prev_pos = snp_positions[i-1]
            
            # Get Transition Matrix based on genetic distance
            trans_matrix = self.transitions.get_transition_matrix(
                genetic_map_func(prev_pos), 
                genetic_map_func(curr_pos)
            )
            trans_matrix = np.maximum(trans_matrix, 1e-300)
            
            curr_em = self.emissions.get_emission_probs(curr_pos, genotypes[i])
            
            for s in range(n_states):
                prev_scores = viterbi[:, i-1] + np.log(trans_matrix[:, s])
                best_prev = int(np.argmax(prev_scores))
                viterbi[s, i] = float(np.max(prev_scores)) + np.log(curr_em[self.states[s]])
                backpointer[s, i] = best_prev

        # 3. Backtrace
        best_path = [np.argmax(viterbi[:, -1])]
        for i in range(n_snps - 1, 0, -1):
            best_path.insert(0, backpointer[best_path[0], i])

        return [self.states[s] for s in best_path]