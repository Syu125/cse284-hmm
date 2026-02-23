import numpy as np

class TransitionModel:
    def __init__(self, generations=10):
        self.G = generations # Generations since admixture

    def get_transition_matrix(self, cm_start, cm_end):
        """
        Calculates the probability of staying in a state vs switching.
        cm_start: genetic position of SNP i
        cm_end: genetic position of SNP i+1
        """
        # 1. Calculate distance in Morgans (cM / 100)
        d = abs(cm_end - cm_start) / 100.0
        
        # 2. Probability of a switch
        # If distance is 0, p_switch is 0.
        p_switch = 1 - np.exp(-self.G * d)
        
        # 3. Probability of staying
        p_stay = 1 - p_switch
        
        # Matrix: [ [Stay, Switch], [Switch, Stay] ]
        # Rows = Current State (YRI, CEU), Cols = Next State (YRI, CEU)
        matrix = np.array([
            [p_stay, p_switch], # From YRI to [YRI, CEU]
            [p_switch, p_stay]  # From CEU to [YRI, CEU]
        ])
        
        return matrix