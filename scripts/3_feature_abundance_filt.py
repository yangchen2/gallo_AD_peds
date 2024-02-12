import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import qiime2 as q2
import logging

# Setup logging
logging.basicConfig(filename='../logs/feature_abundance_filt.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def feature_abundance_filt(biom_path: str, min_samp_prev: float):
    """
    Filter features with at least min_samp_prev prevalence across samples.

    Parameters:
    biom_path (str): The file path to the biom table.
    min_samp_prev (float): Filter features to those with at least X% prevalence across samples

    Returns:
    pd.DataFrame: The dataframe corresponding to the biom table.
    """

    # Load in BIOM file and convert to pandas df
    biom_table = load_table(biom_path)
    df = pd.DataFrame(biom_table.to_dataframe().transpose())
    
    # Convert counts to presence/absence (1/0)
    presence_absence_df = df > 0

    # Calculate prevalence for each feature
    prevalence = presence_absence_df.sum(axis=1) / len(df.columns)

    # Filter features with at least X% prevalence across samples
    filtered_df = df[prevalence >= min_samp_prev].transpose()

    logging.info(f"Features with at least {min_samp_prev} prevalence across samples retained, others filtered out")
    return filtered_df


def save_as_biom(df: pd.DataFrame, output_path: str):
    """
    Save a pandas DataFrame as a BIOM table.

    Parameters:
    df (pd.DataFrame): The DataFrame to save.
    output_path (str): Path to the output BIOM file.
    """
    table = biom.table.Table(df.values, observation_ids=df.index, sample_ids=df.columns)
    with biom_open(output_path, 'w') as f:
        table.to_hdf5(f, "feature abundance filtering")
    logging.info(f"DataFrame saved as BIOM table to {output_path}")


if __name__ == '__main__':
    try:
        biom_path = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare.biom'
        output_path = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare_ab-filt.biom'
        
        # Filter features
        save_as_biom(feature_abundance_filt(biom_path, 0.1), output_path)
        logging.info("Script completed successfully.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
