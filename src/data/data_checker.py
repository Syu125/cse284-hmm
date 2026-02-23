from data_parser import get_population_dict

def main():
    
    # Check sample panel
    panel_path = "../../data/integrated_call_samples_v3.20130502.ALL.panel"
    populations = get_population_dict(panel_path)
    
    print(populations.keys())
    
    print(f"YRI samples: {len(populations.get('YRI', []))}")
    print(f"CEU samples: {len(populations.get('CEU', []))}")
    
if __name__ == "__main__":
    main()