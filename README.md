# Demonstrating the phylotaR pipeline
A series of scripts for demonstrating the phylotaR pipeline for primates and palms. The scripts will first run the phylotaR pipeline and then generate alignments and trees.

## Requirements
- R (version 3+)
- R packages: `phylotaR`, `ape` and `treeman`
- BLAST+
- RAxML
- MAFFT

## Folder set-up
```
- palms
- - [contains all palms phylotaR output]
- primates
- - [contains all primates phylotaR output]
- expected
- - [expected trees of palms and primates]
- figures
- - [tree viz output]
- taxdump.tar.gz
- ncbi-blast-2.7.1+
- - bin
- - - [contains BLAST+ executables]
```

`taxdump.tar.gz` can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
`ncbi-blast-2.7.1+` can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. Download `*-src.tar.gz`, and extract. Note, the version number may be different for you and will need changed in the R scripts. Alternatively these arguments can be left blank in the `run_*.R` scripts and BLAST tools from the system path will be used and phylotaR will download the taxonomy itself.

## Scripts
- `0_run_*.R`: Intiate the phylotaR pipeline.
- `1_cluster_selection_*.R`: Select clusters for alignments and run alignments with MAFFT
- `2_supermatrix_*.R`: Generate supermatrices from alignments and run RAxML
- `3_visualise_*.R`: Visualise the best RAxML tree and compare to expected trees

Each R script must be run manually.

## Figures

Tree generated from two 'best' phylotaR clusters for palms and primates.

#### Palms

![Palms](https://raw.githubusercontent.com/AntonelliLab/phylotaR_demo/master/figures/palms.png)


### Family-level tree comparisons

#### Palms

> Expected tree after Baker et al. (2009) Complete Generic-Level Phylogenetic Analyses of Palms (Arecaceae) with Comparisons of Supertree and Supermatrix Approaches. *Systematic Biology*, 58(2):240–256 [DOI](https://doi.org/10.1093/sysbio/syp021)

![Palms](https://raw.githubusercontent.com/AntonelliLab/phylotaR_demo/master/figures/palms_coplot.png)

