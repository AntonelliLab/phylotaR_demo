# Identify best clusters

# LIBS
library(phylotaR)
source(file.path('tools', 'selection_tools.R'))

# VARS
wd <- 'primates'

# INPUT
clstrs_obj <- genClstrsObj(wd)

# CLTYPE STATS
# drop clusters of 1
n_taxa <- getClNTx(clstrs_obj,
                   clstrs_obj@clstr_ids)
keep <- names(n_taxa)[n_taxa > 10]
clstrs_obj <- drpClstrs(clstrs_obj, keep)
table(sapply(clstrs_obj@clstrs, function(x) x[['cl_type']]))

# FILTER
# get n taxa per cluster
n_taxa <- getClNTx(clstrs_obj,
                   clstrs_obj@clstr_ids)
# drop all clusters with fewer than 100 taxa
keep <- names(n_taxa)[n_taxa > 100]
clstrs_obj <- drpClstrs(clstrs_obj, keep)
for(id in clstrs_obj@clstr_ids) {
  clstrs_obj<- fltrClstrSqs(clstrs_obj, id=id,
                            rank='genus', mn_pambg=0.5)
}
# drop all clusters with MAD scores less than .5
mad_scrs <- getClMAD(clstrs_obj, clstrs_obj@clstr_ids)
keep <- names(mad_scrs)[mad_scrs > .5]
clstrs_obj <- drpClstrs(clstrs_obj, keep)

# SUMMARY
smmry <- genSumTable(clstrs_obj)
write.csv(smmry, file.path('figures', 'best_clusters_primates.csv'))

# SELECT
n_taxa <- getClNTx(clstrs_obj,
                   clstrs_obj@clstr_ids)
keep <- names(n_taxa)[n_taxa > 50]
clstrs_obj <- drpClstrs(clstrs_obj, keep)

# OUTPUT
# write out both clusters
for(i in 1:length(clstrs_obj@clstr_ids)) {
  cid <- clstrs_obj@clstr_ids[[i]]
  getSqDfs(clstrs_obj, id=cid, prse=0.1)
  # get its sequences IDs
  sids <- clstrs_obj@clstrs[[cid]][['gis']]
  # get genus names for fasta def. lines
  scnms <- getIDFrmTxdct(clstrs_obj@txdct, ret='ScientificName',
                         id=names(sids), rank='genus')
  infile <- file.path(wd, paste0('sequences', i, '.fasta'))
  writeSqs(sids=sids, dflns=scnms, flpth=infile)
}

# ALIGN
for(i in 1:length(clstrs_obj@clstr_ids)) {
  inpt <- file.path(wd, paste0('sequences', i,
                               '.fasta'))
  otpt <- file.path(wd, paste0('alignment', i,'.fasta'))
  system(paste0('mafft --auto ', inpt, ' > ', otpt))
}
