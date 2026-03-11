import matplotlib.pyplot as plt

def plot_ancestry(positions, states, true_switch_idx=None, save_path="ancestry_plot.png"):
    """
    Creates a karyogram-style plot of the HMM results.
    """
    # Support both legacy 3-state labels and phased 4-state labels.
    if any(s in {"CEU_YRI", "YRI_CEU"} for s in states):
        val_map = {"YRI_YRI": 0, "YRI_CEU": 1, "CEU_YRI": 2, "CEU_CEU": 3}
        ytick_positions = [0, 1, 2, 3]
        ytick_labels = ["YRI/YRI", "YRI|CEU", "CEU|YRI", "CEU/CEU"]
        y_min, y_max = -0.2, 3.2
    else:
        val_map = {"YRI": 0, "HET": 1, "CEU": 2}
        ytick_positions = [0, 1, 2]
        ytick_labels = ["YRI (African)", "HET (Mixed)", "CEU (European)"]
        y_min, y_max = -0.2, 2.2

    y_values = [val_map[s] for s in states]
    
    plt.figure(figsize=(15, 4))
    
    # Plot the inferred path
    plt.step(positions, y_values, where='post', color='navy', linewidth=2, label='Inferred Ancestry')
    plt.fill_between(positions, y_values, step="post", color='skyblue', alpha=0.3)
    
    # Mark the true switch point if provided
    if true_switch_idx is not None:
        plt.axvline(x=positions[true_switch_idx], color='red', linestyle='--', alpha=0.7, label='Truth Midpoint')

    plt.yticks(ytick_positions, ytick_labels)
    plt.ylim(y_min, y_max)
    plt.title("Local Ancestry Inference: Chromosome 22 (Simulation)")
    plt.xlabel("Physical Position on Chr 22 (bp)")
    plt.legend(loc='center right')
    
    plt.grid(axis='x', linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.savefig(save_path)
    print(f"[+] Plot saved to {save_path}")