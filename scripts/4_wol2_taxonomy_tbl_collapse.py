import pandas as pd
import biom
from biom import load_table
from biom.util import biom_open
import qiime2 as q2
import logging

# Setup logging
logging.basicConfig(filename='../logs/wol2_taxonomy_tbl_collapse.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def group_by_taxonomy_levels(biom_path: str, wol2_taxonomy_path: str, input_filter_type: str):
    """
    Read a rarefied BIOM table, attach gg2 taxonomy info, and collapse by taxonomy level.

    Parameters:
    biom_path (str): Path to the BIOM table
    wol2_taxonomy (str): Path to wol2 taxonomy reference file
    input_filter_type (str): A toggle for whether the input is the abundance filtered table or the zebra filtered table

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
    # Load in wol2 taxonomy file
    wol_taxonomy = pd.read_csv(wol2_taxonomy_path, sep='\t', header=None).set_index(0)
    wol_taxonomy.index.name = None
    wol_taxonomy.rename(columns={wol_taxonomy.columns[0]: 'Taxonomy'}, inplace=True)

    taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    for level in taxonomy_levels:
            # Split Taxon column based on ';'
            wol_taxonomy[taxonomy_levels] = wol_taxonomy['Taxonomy'].str.split(';', expand=True)
    logging.info(f"wol2 taxonomy loaded")

    # Merge taxonomy column based on common index
    df = df.merge(wol_taxonomy, how='left', left_index=True, right_index=True)

    # Define all taxonomy levels
    taxonomy_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    logging.info("STEP 3: Dataframes collapsed by taxonomy level.")
    # Loop through each taxonomy level and perform grouping
    for level in taxonomy_levels:
        # Split Taxon column based on ';'
        df[taxonomy_levels] = df['Taxonomy'].str.split(';', expand=True)

        # Group by (collapse) features to the current taxonomy level
        df_grouped = df.groupby(level).sum(numeric_only=True)

        # Store the grouped dataframe in the dictionary
        grouped_dfs[level] = df_grouped

        # Convert and save dataframes as BIOM tables
        dataframes = {'feature_tables': grouped_dfs[level]}
        for name, table in dataframes.items():
            obs_ids = table.index
            samp_ids = table.columns

            biom_table = biom.table.Table(table.values, observation_ids=obs_ids, sample_ids=samp_ids)
            biom_output_file = f"../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare_{input_filter_type}-filt_{level}_collapsed.biom"

            with biom_open(biom_output_file, 'w') as f:
                biom_table.to_hdf5(f, generated_by="collapsed tables by taxonomy level")  
    
    logging.info("STEP 4: Tables collapsed at various taxonomy levels.")
    return grouped_dfs


if __name__ == '__main__':
    try:
        # File paths for the BIOM files
        biom_path_abundance = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare_ab-filt.biom'
        biom_path_zebra = '../tables/fastp_hg38_t2t_pangenome_193238_feature-table_rare_zebra-filt.biom'
        wol2_taxonomy_path = '../wol2_refs/lineages.txt'

        # Create collapsed taxa BIOMs
        logging.info("-----> Performing script with abundance filtered input.")
        group_by_taxonomy_levels(biom_path_abundance, wol2_taxonomy_path, "ab")

        logging.info("-----> Performing script with Zebra filtered input.")
        group_by_taxonomy_levels(biom_path_zebra, wol2_taxonomy_path, "zebra")

        logging.info("Script completed successfully.")

    except Exception as e:
        logging.error(f"An error occurred in the main execution: {e}")