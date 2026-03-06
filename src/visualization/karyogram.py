import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

def plot_ancestry(positions, states, true_switch_idx=None, save_path="ancestry_plot.png"):
    """
    Creates a karyogram-style plot of the HMM results.
    """
    val_map = {"YRI": 0, "HET": 1, "CEU": 2}
    y_values = [val_map[s] for s in states]
    
    plt.figure(figsize=(15, 4))
    
    # Plot the inferred path
    plt.step(positions, y_values, where='post', color='navy', linewidth=2, label='Inferred Ancestry')
    plt.fill_between(positions, y_values, step="post", color='skyblue', alpha=0.3)
    
    # Mark the true switch point if provided
    if true_switch_idx is not None:
        plt.axvline(x=positions[true_switch_idx], color='red', linestyle='--', alpha=0.7, label='Truth Midpoint')

    plt.yticks([0, 1, 2], ["YRI (African)", "HET (Mixed)", "CEU (European)"])
    plt.ylim(-0.2, 2.2)
    plt.title("Local Ancestry Inference: Chromosome 22 (Simulation)")
    plt.xlabel("Physical Position on Chr 22 (bp)")
    plt.legend(loc='center right')
    
    plt.grid(axis='x', linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.savefig(save_path)
    print(f"[+] Plot saved to {save_path}")