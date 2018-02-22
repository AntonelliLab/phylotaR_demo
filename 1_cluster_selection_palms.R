# Identify best clusters

# LIBS
library(phylotaR)
source(file.path('tools', 'selection_tools.R'))
source('dev.R')

# VARS
wd <- 'palms'

# INPUT
clstrs_obj <- genClstrsObj(wd)

# CLTYPE STATS
# drop clusters of 1
n_taxa <- getClNTx(clstrs_obj,
                   clstrs_obj@clstr_ids)
keep <- names(n_taxa)[n_taxa > 10]
clstrs_obj <- drpClstrs(clstrs_obj, keep)
table(sapply(clstrs_obj@clstrs, function(x) x[['cl_type']]))

# REDUCE
# get n taxa per cluster
n_taxa <- getClNTx(clstrs_obj,
                   clstrs_obj@clstr_ids)
# drop all clusters with fewer than 100 taxa
keep <- names(n_taxa)[n_taxa > 100]
clstrs_obj <- drpClstrs(clstrs_obj, keep)

# DROP SUSPECT SEQUENCES
blcklst <- c(62511848, 353531054, 440579093)
clstrs_obj <- rmSqs(clstrs_obj, blcklst)

# FILTER
for(id in clstrs_obj@clstr_ids) {
  clstrs_obj<- fltrClstrSqs(clstrs_obj, id=id,
                            rnk='tribe', mn_pambg=0.5,
                            n_ech=2)
}

# REDUCE
# drop all clusters with MAD scores less than .5
mad_scrs <- getClMAD(clstrs_obj, clstrs_obj@clstr_ids)
keep <- names(mad_scrs)[mad_scrs > .5]
clstrs_obj <- drpClstrs(clstrs_obj, keep)

# SUMMARY
smmry <- genSumTable(clstrs_obj)
write.csv(smmry, file.path('figures', 'best_clusters_palms.csv'))

# OUTPUT
# write out top 10 clusters with most taxa
cids <- smmry[['id']][1:10]
sqfls <- list.files(wd, pattern='sequences[0-9]+.fasta')
rmFls(file.path(wd, sqfls))
for(i in seq_along(cids)) {
  cid <- clstrs_obj@clstr_ids[[i]]
  getSqDfs(clstrs_obj, id=cid, prse=0.1)
  # get its sequences IDs
  sids <- clstrs_obj@clstrs[[cid]][['gis']]
  tids <- sapply(sids, function(x) clstrs_obj@sqs[[x]][['ti']])
  # get sci names for fasta def. lines
  scnms <- getIDFrmTxdct(clstrs_obj@txdct, ret='ScientificName',
                         id=tids, rnk='tribe')
  n <- sapply(seq_along(scnms), function(x) 
    sum(scnms[x] == scnms[x:length(scnms)]))
  infile <- file.path(wd, paste0('sequences', i, '.fasta'))
  writeSqs(sids=sids, dflns=paste0(scnms, '_', n), flpth=infile)
}

# ALIGN
alfls <- list.files(wd, pattern='alignment[0-9]+.fasta')
rmFls(file.path(wd, alfls))
for(i in seq_along(cids)) {
  inpt <- file.path(wd, paste0('sequences', i,
                               '.fasta'))
  otpt <- file.path(wd, paste0('alignment', i,'.fasta'))
  system(paste0('mafft --auto ', inpt, ' > ', otpt))
}
