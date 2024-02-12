import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import numpy as np
from numpy.random import RandomState
import logging

# Setup logging
logging.basicConfig(filename='../logs/rarefy.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def read_and_convert_biom(biom_path: str):
    """
    Read a BIOM-format table and convert it into a Pandas dataframe.

    Parameters:
    biom_path (str): The file path to the biom table.

    Returns:
    pd.DataFrame: The dataframe corresponding to the biom table.
    """
    logging.info(f"Loading BIOM file: {biom_path}")
    try:
        biom_table = load_table(biom_path)
        df = pd.DataFrame(biom_table.to_dataframe().T)
        logging.info(f"BIOM table shape: {df.shape}")
        logging.info("BIOM file successfully converted to DataFrame")
        return df
    except Exception as e:
        logging.error(f"Error in processing BIOM file: {e}")
        raise

def rarefy_table(biom_df: pd.DataFrame, seed: int = 42) -> pd.DataFrame:
    """
    Rarefy a single BIOM table to its minimum sampling depth.

    Parameters:
    biom_df (pd.DataFrame): Input DataFrame representing the BIOM table.
    seed (int): Seed for the random number generator to ensure reproducibility.

    Returns:
    pd.DataFrame: A rarefied DataFrame.
    """
    # Calculate minimum sampling depth
    # min_depth = biom_df.sum(axis=1).min().astype(int)
    min_depth = 1000 # Set min depth to 1000
    logging.info(f"Minimum sampling depth: {min_depth}")

    # Rarefaction process
    prng = RandomState(seed)
    rarefied_data = np.zeros_like(biom_df.values)

    for i, row in enumerate(biom_df.values):
        total_count = row.sum()
        if total_count >= min_depth:
            probabilities = row / total_count
            chosen_indices = prng.choice(biom_df.columns.size, min_depth, p=probabilities)
            rarefied_row = np.bincount(chosen_indices, minlength=biom_df.columns.size)
            rarefied_data[i] = rarefied_row
    
    rarefied_df = pd.DataFrame(rarefied_data, index=biom_df.index, columns=biom_df.columns)
    logging.info("Rarefaction completed.")
    return rarefied_df

def save_as_biom(df: pd.DataFrame, output_path: str):
    """
    Save a pandas DataFrame as a BIOM table.

    Parameters:
    df (pd.DataFrame): The DataFrame to save.
    output_path (str): Path to the output BIOM file.
    """
    table = biom.table.Table(df.values, observation_ids=df.index, sample_ids=df.columns)
    with biom_open(output_path, 'w') as f:
        table.to_hdf5(f, "rarefaction script")
    logging.info(f"DataFrame saved as BIOM table to {output_path}")

if __name__ == '__main__':
    try:
        biom_path = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table.biom'
        output_path = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare.biom'
        
        df = read_and_convert_biom(biom_path)
        rarefied_df = rarefy_table(df)
        save_as_biom(rarefied_df, output_path)
        
        logging.info("Script completed successfully.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
