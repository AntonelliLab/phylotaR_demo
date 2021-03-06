# Demonstrating the phylotaR pipeline
A series of scripts for demonstrating the phylotaR pipeline for primates and palms. The scripts will first run the phylotaR pipeline and then generate alignments and trees.

## Requirements
- R (version 3+)
- R packages: `phylotaR`, `ape` and `treeman`
- BLAST+
- RAxML
- MAFFT

## Install specific phylotaR release

```{r}
# this will install the exact development version used to generate results
devtools::install_github(repo='AntonelliLab/phylotaR', ref='5677e1560f7b0f8f60e5109072d5af0326338d69')
# the code is not guarranteed to work with any other version
```

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

Genus-level trees generated from the two 'best' phylotaR clusters for palms and primates.

### Palms

![Palms](https://raw.githubusercontent.com/AntonelliLab/phylotaR_demo/master/figures/palms.png)

**Family-level comparison**
> Expected tree from Baker et al. (2009) Complete Generic-Level Phylogenetic Analyses of Palms (Arecaceae) with Comparisons of Supertree and Supermatrix Approaches. *Systematic Biology*, 58(2):240–256 [DOI](https://doi.org/10.1093/sysbio/syp021)

![Palms](https://raw.githubusercontent.com/AntonelliLab/phylotaR_demo/master/figures/palms_coplot.png)

### Primates

![Primates](https://raw.githubusercontent.com/AntonelliLab/phylotaR_demo/master/figures/primates.png)

**Comparison**
> Expected tree from Perelman et al. (2011) A Molecular Phylogeny of Living Primates. *PLOS Genetics* [DOI](https://doi.org/10.1371/journal.pgen.1001342)

![Primates](https://raw.githubusercontent.com/AntonelliLab/phylotaR_demo/master/figures/primates_coplot.png)

