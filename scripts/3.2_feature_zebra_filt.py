import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import qiime2 as q2
import logging

# Setup logging
logging.basicConfig(filename='../logs/feature_zebra_filt.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def feature_abundance_filt(biom_path: str, zebra_coverage_path: str, zebra_threshold: float, filtered_biom_path: str):
    """
    Filter out features  below a certain Zebra filter threshold.

    Parameters:
    biom_path (str): The file path to the biom table
    zebra_coverage_path (str): The file path to the zebra coverages
    zebra_threshold (float): Float while where OTUs below that threshold will be filtered out
    filtered_biom_path (str): The file path of the outputted zebra filtered biom table

    Returns:
    pd.DataFrame: The dataframe corresponding to the zebra filterd biom table.
    """

    # Load in BIOM file and convert to pandas df
    logging.info("Step 1: Loading in biom table and converting to pandas df")
    biom_table = load_table(biom_path)
    df = pd.DataFrame(biom_table.to_dataframe().transpose())
    logging.info(f"Table shape: {df.shape}")

    # Load in zebra coverages file
    logging.info("Step 2: Loading in Zebra filter coverages output file")
    zebra_coverage = pd.read_csv(zebra_coverage_path, sep='\t', header=None).set_index(0)
    zebra_coverage.index.name = None
    zebra_coverage.rename(columns={zebra_coverage.columns[0]: 'coverage %'}, inplace=True)
    zebra_coverage_sorted = zebra_coverage.sort_values(by='coverage %', ascending=False)
    logging.info("Saving Zebra filter coverages file as sorted")
    zebra_coverage_sorted.to_csv('../zebra_coverage/zebra_coverage_sorted.tsv') 

    logging.info("Step 3: Filtering df based on the coverage percentage condition")
    # Create a boolean mask where True represents rows with coverage >= 10.00%
    coverage_mask = zebra_coverage_sorted['coverage %'] >= zebra_threshold
    # Use the mask to filter 'zebra_coverage_sorted' first to get the list of Feature OTUs meeting the condition
    filtered_features = zebra_coverage_sorted[coverage_mask].index
    # Now, filter 'df' to keep only the rows where the Feature OTUs are in 'filtered_features'  
    filtered_df = df[df.index.isin(filtered_features)]
    logging.info(f"Table shape after zebra filtering: {filtered_df.shape}")
    filtered_df = filtered_df.transpose()
    
    logging.info("Step 4: Saving filtered table as a new BIOM file")
    table = biom.table.Table(filtered_df.values, observation_ids=filtered_df.index, sample_ids=filtered_df.columns)
    with biom_open(filtered_biom_path, 'w') as f:
        table.to_hdf5(f, "feature abundance filtering")
    logging.info(f"DataFrame saved as BIOM table to {filtered_biom_path}")



if __name__ == '__main__':
    try:
        biom_path = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare.biom'
        zebra_coverage_path = '../zebra_coverage/coverage_percentage.txt'
        filtered_biom_path = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare_zebra-filt.biom'
        
        # Filter features
        feature_abundance_filt(biom_path, zebra_coverage_path, 10.0, filtered_biom_path)
        logging.info("Script completed successfully.")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise