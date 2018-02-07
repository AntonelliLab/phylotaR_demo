# Demonstrating the phylotaR pipeline
A series of scripts for demonstrating the phylotaR pipeline for primates and palms.

## Requirements
- R (version 3+)
- phylotaR package
- BLAST+
- RAxML
- MAFFT

## Scripts
- `run_*.R`: Intiate the phylotaR pipeline.
- `cluster_selection_*.R`: Select clusters for alignments
- `align_*.sh`: Align selected clusters with MAFFT
- `supermatrix_*.R`: Generate supermatrices from alignments
- `construct_*.sh`: Construct tree from supermatrix with RAxML
- `visualise_*.R`: Visualise the best RAxML tree

Each R script must be run manually. Shell scripts are called from within R.
