import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import qiime2 as q2
import logging

# Setup logging
logging.basicConfig(filename='../logs/gg_taxonomy_tbl_collapse.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Load in Greengenes2 taxonomy file
logging.info(f"Loading GreenGreens2 taxonomy... (this will take 1-2 minutes)")
gg_taxonomy = q2.Artifact.load('/Users/yac027/gg_taxonomy/2022.10.taxonomy.asv.tsv.qza').view(pd.DataFrame)
logging.info(f"GreenGreens2 taxonomy loaded")


def group_by_taxonomy_levels(biom_path: str, gg_taxonomy: pd.DataFrame):
    """
    Read a rarefied BIOM table, attach gg2 taxonomy info, and collapse by taxonomy level.

    Parameters:
    biom_path (str): Path to the BIOM table
    gg_taxonomy (pd.DataFrame): DataFrame with taxonomy information from Gg2

    Returns:
    dict: A dictionary containing dataframes collapsed at each of 7 taxonomy levels.
    """

    logging.info("STEP 1: Converting BIOM to dataframe.")
    # Read in BIOM table and convert to Pandas df
    biom_table = load_table(biom_path)
    df = pd.DataFrame(biom_table.to_dataframe()).transpose()

    # Dictionary to hold the grouped dataframes for each level
    grouped_dfs = {}

    logging.info("STEP 2: Attaching taxonomy info and collapsing by taxonomy level.")
    # Merge taxonomy column based on common index
    df = df.merge(gg_taxonomy, how='left', left_index=True, right_index=True)

    # Define all taxonomy levels
    taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    logging.info("STEP 3: Dataframes collapsed by taxonomy level.")
    # Loop through each taxonomy level and perform grouping
    for level in taxonomy_levels:
        # Split Taxon column based on ';'
        df[taxonomy_levels] = df['Taxon'].str.split(';', expand=True)

        # Group by (collapse) features to the current taxonomy level
        df_grouped = df.groupby(level).sum()

        # Store the grouped dataframe in the dictionary
        grouped_dfs[level] = df_grouped

    logging.info("STEP 4: Tables collapsed at various taxonomy levels.")
    # Convert and save dataframes as BIOM tables
    dataframes = {'feature_tables': grouped_dfs[level]}
    for name, table in dataframes.items():
        obs_ids = table.index
        samp_ids = table.columns

        biom_table = biom.table.Table(table.values, observation_ids=obs_ids, sample_ids=samp_ids)
        biom_output_file = f"../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare_ab-filt_{level}_collapsed.biom"

        with biom_open(biom_output_file, 'w') as f:
            biom_table.to_hdf5(f, generated_by="collapsed tables by taxonomy level")  
    
    return grouped_dfs


if __name__ == '__main__':
    try:
        # File paths for the BIOM files
        biom_path = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare_ab-filt.biom'

        # Create collapsed taxa BIOMs
        group_by_taxonomy_levels(biom_path, gg_taxonomy)
        logging.info("Script completed successfully.")

    except Exception as e:
        logging.error(f"An error occurred in the main execution: {e}")