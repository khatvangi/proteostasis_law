# Save this as download_genome.py and run with 'python download_genome.py'
from datasets import load_dataset_builder
import os

accession_id = "GCF_000005845.2"
output_dir = "ecoli_MG1655"
output_zip = f"{output_dir}.zip"

# NOTE: The 'datasets' library is for ML data, not general genome downloads.
# The original command 'datasets download genome accession...' suggests
# you might be using a specific, niche tool that happens to be named 'datasets'.

# If you are using the NCBI Datasets CLI tool specifically:
# The correct tool is usually installed via Conda as 'ncbi-datasets-cli'.

# Try installing the correct NCBI tool if the 'datasets' in your command refers to NCBI
# conda install -c bioconda ncbi-datasets-cli

