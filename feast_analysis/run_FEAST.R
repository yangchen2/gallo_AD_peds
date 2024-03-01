# This script runs FEAST: https://github.com/cozygene/FEAST

setwd("/Users/yac027/Gallo_lab/gallo_AD_peds/feast_analysis/")

# Load inputs
metadata <- Load_metadata(metadata_path = "data/Streptococcus_metadata.tsv")
otus <- Load_CountMatrix(CountMatrix_path = "data/Streptococcus_otu.tsv")

# Run FEAST
FEAST_output <- FEAST(C = otus, metadata = metadata, 
                      different_sources_flag = 0, 
                      dir_path = "data",
                      outfile="strep_oral-source")
