import numpy as np

class InferenceEngine:
    def __init__(self, emission_model, transition_model):
        self.emissions = emission_model
        self.transitions = transition_model
        self.states = ["YRI", "CEU"]

    def run_viterbi(self, snp_positions, genotypes, genetic_map_func):
        """
        snp_positions: list of physical positions
        genotypes: list of tuples (0,1), (1,1), etc.
        genetic_map_func: a lambda or function that takes pos and returns cM
        """
        n_snps = len(snp_positions)
        n_states = len(self.states)
        
        # viterbi[state, snp] stores the max probability to reach that point
        viterbi = np.zeros((n_states, n_snps))
        # backpointer stores which state we came from
        backpointer = np.zeros((n_states, n_snps), dtype=int)

        # 1. Initialize (First SNP)
        first_pos = snp_positions[0]
        first_em = self.emissions.get_emission_probs(first_pos, genotypes[0])
        viterbi[0, 0] = np.log(0.5 * first_em["YRI"]) # Log math prevents underflow
        viterbi[1, 0] = np.log(0.5 * first_em["CEU"])

        # 2. Recursion
        for i in range(1, n_snps):
            curr_pos = snp_positions[i]
            prev_pos = snp_positions[i-1]
            
            # Get Transition Matrix based on genetic distance
            trans_matrix = self.transitions.get_transition_matrix(
                genetic_map_func(prev_pos), 
                genetic_map_func(curr_pos)
            )
            
            curr_em = self.emissions.get_emission_probs(curr_pos, genotypes[i])
            
            for s in range(n_states):
                # Calculate: Prev_Viterbi + Log_Transition + Log_Emission
                # We use Log space because multiplying 10,000 small probabilities 
                # would result in a number so small a computer can't store it.
                prob_from_yri = viterbi[0, i-1] + np.log(trans_matrix[0, s])
                prob_from_ceu = viterbi[1, i-1] + np.log(trans_matrix[1, s])
                
                # Pick the best path to this state
                best_prev = np.argmax([prob_from_yri, prob_from_ceu])
                viterbi[s, i] = max(prob_from_yri, prob_from_ceu) + np.log(curr_em[self.states[s]])
                backpointer[s, i] = best_prev

        # 3. Backtrace
        best_path = [np.argmax(viterbi[:, -1])]
        for i in range(n_snps - 1, 0, -1):
            best_path.insert(0, backpointer[best_path[0], i])

        return [self.states[s] for s in best_path]