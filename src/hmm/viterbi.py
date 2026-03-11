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
        self.states = list(getattr(emission_model, "states", ["CEU", "YRI"]))

    def run_viterbi(self, snp_positions, alleles, genetic_map_func):
        """
        snp_positions: list of physical positions
        alleles: list of haplotype alleles (0 or 1)
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
        first_em = self.emissions.get_haplotype_emission_probs(first_pos, alleles[0])
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
            
            curr_em = self.emissions.get_haplotype_emission_probs(curr_pos, alleles[i])
            
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

    def run_viterbi_phased(self, snp_positions, genotypes, genetic_map_func):
        """
        Runs haplotype-level Viterbi on both phased alleles and returns two paths.
        genotypes: list of tuples (a1, a2)
        """
        hap1 = [gt[0] for gt in genotypes]
        hap2 = [gt[1] for gt in genotypes]
        hap1_states = self.run_viterbi(snp_positions, hap1, genetic_map_func)
        hap2_states = self.run_viterbi(snp_positions, hap2, genetic_map_func)
        return hap1_states, hap2_states

    @staticmethod
    def combine_haplotype_states(hap1_states, hap2_states):
        return [f"{h1}_{h2}" for h1, h2 in zip(hap1_states, hap2_states)]