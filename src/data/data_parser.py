"""
This module contains functions for parsing genetic data.
"""

import pandas as pd

def get_population_dict(panel_file):
    df = pd.read_csv(panel_file, sep='\t')
    populations = df.groupby('pop')['sample'].apply(list).to_dict()
    return populations

