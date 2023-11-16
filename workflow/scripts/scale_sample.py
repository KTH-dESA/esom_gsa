"""Module for scaling a morris sample"""
import sys 
import numpy as np
from typing import Dict
import pandas as pd

def main(morris_sample: np.array, parameters: pd.DataFrame, value_for_interpolate: str = "end") -> np.array:
    """Very important to preserve order of sample and parameters"""
    
    scaled = morris_sample.copy()
    shape = morris_sample.shape
    
    if value_for_interpolate == "start": # use base year values
        for col in range(shape[1]):
            difference = parameters.loc[col, "max_value_base_year"] - parameters.loc[col, "min_value_base_year"]
            scaled[:,col] = scaled[:,col] * difference + parameters.loc[col, "min_value_end_year"]
    elif value_for_interpolate == "end": # use end year values
        for col in range(shape[1]):
            difference = parameters.loc[col, "max_value_end_year"] - parameters.loc[col, "min_value_end_year"]
            scaled[:,col] = scaled[:,col] * difference + parameters.loc[col, "min_value_end_year"]
    else: # use average 
        for col in range(shape[1]):
            scaled_base = morris_sample.copy()
            difference_base = parameters.loc[col, "max_value_base_year"] - parameters.loc[col, "min_value_base_year"]
            scaled_base[:,col] = scaled_base[:,col] * difference_base + parameters.loc[col, "min_value_base_year"]
            
            scaled_end = morris_sample.copy()
            difference_end = parameters.loc[col, "max_value_end_year"] - parameters.loc[col, "min_value_end_year"]
            scaled_end[:,col] = scaled_end[:,col] * difference_end + parameters.loc[col, "min_value_end_year"]
            
            scaled[:,col] = (scaled_base[:,col] + scaled_end[:,col]) / 2
    
    return scaled
    
if __name__ == "__main__":
    
    sample_file = sys.argv[1]
    parameters_file = sys.argv[2]
    scaled_file = sys.argv[3]
    params = pd.read_csv(parameters_file)
    morris_sample = np.loadtxt(sample_file, delimiter=",")
    scaled = main(morris_sample, params)
    np.savetxt(scaled_file, scaled, delimiter=',')