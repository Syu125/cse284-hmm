import numpy as np

class EmissionModel:
    def __init__(self, yri_freqs, ceu_freqs):
        """
        yri_freqs: dict {pos: freq}
        ceu_freqs: dict {pos: freq}
        """
        self.yri_freqs = yri_freqs
        self.ceu_freqs = ceu_freqs

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
        return prob