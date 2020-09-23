"""Example data processing script
"""
import sys

import pandas as pd

def read_data(list_of_filepaths):
    df = pd.read_csv()

    return df

def write_data(df, list_of_filepaths):
    """Write out data to disk
    """
    pass

def process_data(df):
    """
    """
    

    return processed_df


def main(list_of_input_filepaths, list_of_output_filepaths):

    df = read_data(list_of_input_filepaths)

    processed_data = process_data(df)

    write_data(process_data, list_of_output_filepaths)

if __name__ == "__main__":

    main(sys.argv[1], sys.argv[2])