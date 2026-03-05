# Implementation Details

This document describes the technical implementation of the HMM for local ancestry inference.

---

## Project Overview

**Problem**: Determine the ancestry (YRI vs CEU) at each location along the genome for an admixed individual.

**Solution**: Use a Hidden Markov Model where:
- **Hidden states**: Ancestry at each SNP (YRI or CEU)
- **Observations**: Genotypes at each SNP
- **Goal**: Find most likely ancestry sequence given observed genotypes

---

## HMM Architecture

### States
Two possible states at each SNP:
- State 0: YRI ancestry
- State 1: CEU ancestry

### Observations
Genotypes for each SNP: {0, 1, 2} (number of alternate alleles)

### Model Components

```python
├── EmissionModel()
│   ├── Input: genotype, allele frequency, state
│   └── Output: P(genotype | state)
│
├── TransitionModel()
│   ├── Input: genetic distance, generations since admixture
│   └── Output: P(state[i+1] | state[i])
│
└── ViterbiAlgorithm()
    ├── Input: emissions, transitions
    └── Output: Most likely state sequence
```

---

## Emission Model

**Purpose**: Calculate $P(\text{genotype} | \text{ancestry state})$

### Hardy-Weinberg Equilibrium Assumption

For a population with allele frequency $p$:

$$P(\text{genotype} = g | \text{state}) = P(a_1 | \text{state}) \times P(a_2 | \text{state})$$

where $a_1, a_2$ are the two alleles and:
- $P(\text{allele} = 0 | \text{state}) = 1 - p$
- $P(\text{allele} = 1 | \text{state}) = p$

### Genotype Probabilities

| Genotype | Alleles | Probability |
|----------|---------|-------------|
| 0 | 0/0 | $(1-p)^2$ |
| 1 | 0/1 or 1/0 | $2p(1-p)$ |
| 2 | 1/1 | $p^2$ |

### Implementation

```python
# src/hmm/emission.py
def emission_prob(genotype, allele_freq, state):
    """
    Calculate P(genotype | state)
    
    Args:
        genotype: 0, 1, or 2 (number of alternate alleles)
        allele_freq: p (frequency of alternate allele)
        state: 0 (YRI) or 1 (CEU)
    
    Returns:
        log probability
    """
    p = yri_freq if state == 0 else ceu_freq
    
    if genotype == 0:
        return log((1-p)**2)
    elif genotype == 1:
        return log(2*p*(1-p))
    else:  # genotype == 2
        return log(p**2)
```

### Allele Frequency Computation

For each SNP, calculate frequency from reference population:

$$p = \frac{\text{# alternate alleles in population}}{2 \times \text{# individuals}}$$

Example: If 100 YRI individuals, 80 with genotype 1/1 (160 alt alleles), 20 with 0/0:
$$p_{YRI} = \frac{160}{200} = 0.80$$

---

## Transition Model

**Purpose**: Calculate $P(\text{state}[i+1] | \text{state}[i])$ based on recombination

### Recombination Probability

Under simple admixture model with $G$ generations since admixture:

$$P(\text{switch at genetic distance } d) = 1 - e^{-G \cdot d}$$

where:
- $d$ = genetic distance between consecutive SNPs (in Morgans)
- $G$ = generations since admixture (typical: 10-20 for recent admixture)
- Assumes one crossover per meiosis = 1 Morgan

### Transition Matrix

$$T = \begin{pmatrix} 1 - P_{switch} & P_{switch} \\ P_{switch} & 1 - P_{switch} \end{pmatrix}$$

### Implementation

```python
# src/hmm/transition.py
def transition_prob(genetic_distance, generations=10):
    """
    Probability of ancestry switch between consecutive SNPs
    
    Args:
        genetic_distance: Distance in Morgans (cM / 100)
        generations: Time since admixture
    
    Returns:
        Probability of state change [0, 1]
    """
    p_switch = 1 - np.exp(-generations * genetic_distance)
    return p_switch
```

### Genetic Map Interpolation

Convert physical position (bp) to genetic position (cM):

```python
# src/data/data_parser.py
def interpolate_genetic_position(physical_pos, physical_map, genetic_map):
    """
    Linearly interpolate physical position to genetic position
    """
    idx = np.searchsorted(physical_map, physical_pos)
    
    if idx == 0:
        return genetic_map[0]
    elif idx == len(physical_map):
        return genetic_map[-1]
    else:
        # Linear interpolation
        w = (physical_pos - physical_map[idx-1]) / \
            (physical_map[idx] - physical_map[idx-1])
        return genetic_map[idx-1] + w * (genetic_map[idx] - genetic_map[idx-1])
```

---

## Viterbi Algorithm

**Purpose**: Find most likely state sequence given observations

### Problem
Given observations $O = (o_1, o_2, ..., o_n)$ and model parameters, find:

$$\text{argmax}_S \prod_{i=1}^{n} P(o_i | s_i) \cdot P(s_i | s_{i-1})$$

### Solution: Dynamic Programming

Define:
$$\delta_i(s) = \max_{s_1...s_{i-1}} \prod_{j=1}^{i} P(o_j | s_j) \cdot P(s_j | s_{j-1})$$

**Recurrence**:
$$\delta_{i+1}(s') = P(o_{i+1} | s') \cdot \max_s [\delta_i(s) \cdot P(s' | s)]$$

### Implementation (Log-Space)

Work in log-space to avoid numerical underflow:

$$\log \prod x_i = \sum \log x_i$$

```python
# src/hmm/viterbi.py
def viterbi(emissions, transitions):
    """
    Viterbi decoding in log space
    
    Args:
        emissions: (n_snps, n_states) matrix of log P(obs | state)
        transitions: (n_states, n_states) matrix of log P(state' | state)
    
    Returns:
        path: (n_snps,) most likely state sequence
    """
    n_snps, n_states = emissions.shape
    
    # Initialize
    viterbi_probs = np.zeros((n_snps, n_states))
    viterbi_path = np.zeros((n_snps, n_states), dtype=int)
    
    viterbi_probs[0] = emissions[0] + np.log(1 / n_states)  # Uniform prior
    
    # Forward pass
    for i in range(1, n_snps):
        for s in range(n_states):
            # Max over previous states
            prev_max = viterbi_probs[i-1] + transitions[:, s]
            viterbi_path[i, s] = np.argmax(prev_max)
            viterbi_probs[i, s] = emissions[i, s] + np.max(prev_max)
    
    # Backtrack
    path = np.zeros(n_snps, dtype=int)
    path[-1] = np.argmax(viterbi_probs[-1])
    for i in range(n_snps-2, -1, -1):
        path[i] = viterbi_path[i+1, path[i+1]]
    
    return path
```

### Time Complexity
- Training: $O(n \times K^2)$ where $n$ = # SNPs, $K$ = # states
- With K=2 (YRI/CEU): $O(4n)$ = linear
- For 500k SNPs: ~2 million operations

---

## Data Loading Pipeline

### VCF Parsing

Extract genotypes and compute allele frequencies:

```python
# src/data/data_parser.py
def get_allele_frequencies(vcf_path, sample_indices):
    """
    Read VCF file and compute allele frequencies for samples
    
    Returns: dict mapping SNP position to allele frequency
    """
    vcf = pysam.VariantFile(vcf_path)
    freqs = {}
    
    for record in vcf:
        total_alleles = 0
        alt_alleles = 0
        
        for sample_idx in sample_indices:
            gt = record.samples[sample_idx]['GT']
            alt_count = sum(1 for allele in gt if allele == 1)
            alt_alleles += alt_count
            total_alleles += 2  # Diploid
        
        # Frequency = alt_alleles / total_alleles
        freqs[record.pos] = alt_alleles / total_alleles
    
    return freqs
```

### Sample Panel Parsing

Extract population assignments:

```python
def get_population_dict(panel_path):
    """
    Parse 1000 Genomes panel file
    
    Returns: dict mapping population to sample indices
    """
    pop_dict = defaultdict(list)
    
    with open(panel_path) as f:
        next(f)  # Skip header
        for idx, line in enumerate(f):
            parts = line.strip().split()
            pop = parts[1]  # Column 2 is population code
            pop_dict[pop].append(idx)
    
    return pop_dict
```

---

## Analysis Workflow

### 1. Real Sample Analysis (Production)
```
Input: Reference VCF + Query VCF
├─ Extract YRI/CEU reference samples
├─ Compute allele frequencies
├─ For each ASW sample:
│  ├─ Create emission probabilities
│  ├─ Run Viterbi
│  └─ Store results
└─ Output: Results table + plots
```

---

## Performance Considerations

### Memory Usage
```
YRI frequencies:    ~500k SNPs × 8 bytes = 4 MB
CEU frequencies:    ~500k SNPs × 8 bytes = 4 MB
Emission matrix:    ~500k SNPs × 2 states × 8 bytes = 8 MB
Viterbi DP table:   ~500k SNPs × 2 states × 8 bytes = 8 MB
─────────────────────────────────────────────────
Total per sample:   ~24 MB
All 1000+ samples:  ~24 GB (problematic!)
```

**Solutions**:
- Cache frequencies (reuse for all samples)
- Process samples sequentially (one at a time)
- Use slice dataset for memory-constrained machines

### Runtime Analysis
```
Frequency calculation: 5-10 min (first run, full dataset)
Frequency caching:     <1 sec (cached)
Viterbi per sample:    100-500 ms
─────────────────────────────────────────────
Total (100 samples):   ~10-50 minutes
```

---

## Key Assumptions

1. **Hardy-Weinberg Equilibrium**: Random mating within populations
2. **Linkage Independence**: Genotypes only dependent on ancestry state (no LD)
3. **Known Admixture Time**: Generations since admixture is constant
4. **Simple Two-Population Model**: Real admixture may involve >2 sources
5. **No Assortative Mating**: Ancestry doesn't correlate with partner choice

**Real-world violations**:
- LD causes non-independence → affects confidence intervals
- Admixture may be multi-way (YRI + CEU + Native American)
- Variable admixture time across genome

---

## References

1. **Viterbi Algorithm**: Viterbi, A. (1967). Error bounds for convolutional codes. *IEEE Transactions on Information Theory*
2. **Admixture Mapping**: McKeigue, P. M. (1998). Mapping genes for complex diseases. *American Journal of Human Genetics*
3. **Local Ancestry**: Price, A. L., et al. (2010). Sensitive detection of ancestry inference. *PLoS Genetics*
4. **1000 Genomes**: The 1000 Genomes Project Consortium. (2015). *Nature*, 526, 68-74

---

## Extending the Model

### Multiple Populations
Replace binary states with K populations (YRI, CEU, Native American, etc.)
- Transition matrix: $K \times K$ instead of $2 \times 2$
- Emission calculation: Same HWE formula, but for each population

### Variable Admixture Time
Estimate generations from data using EM algorithm
- Currently hardcoded: `generations = 10`
- Can be inferred from switch rate distribution

### Population-Specific Recombination Rates
Use population-specific genetic maps
- Currently: Common genetic map for all populations
- Could use YRI, CEU, or other population-specific maps

### Confidence Intervals
Use posterior probabilities instead of Viterbi:
- Compute all possible paths with likelihoods
- Posterior: P(state | data) for each state
- Confidence: Posterior probability of most likely state

