# Identify best clusters

# LIBS
library(phylotaR)
source(file.path('tools', 'selection_tools.R'))

# VARS
wd <- 'primates'

# INPUT
all_cls <- read_phylota(wd)

# CLTYPE STATS
# drop clusters of 10
# n_taxa <- get_ntaxa(all_cls, cid=all_cls@cids)
# n_taxa <- get_cl_slot(all_cls, cid=all_cls@cids, slt_nm='ntx')
# keep <- names(n_taxa)[n_taxa > 10]
keep <- all_cls@cids[1:100]
all_cls <- drop_cls(all_cls, keep)
table(sapply(all_cls@cls@cls, function(x) x@typ))

# REDUCE
# get n taxa per cluster
#n_taxa <- get_ntaxa(all_cls, cid=all_cls@cids)
# drop all clusters with fewer than 100 taxa
#keep <- names(n_taxa)[n_taxa > 100]
#all_cls <- drop_cls(all_cls, keep)

# FILTER
genus_only <- drop_by_rank(all_cls, rnk='genus', n=2,
                           choose_by=c('nncltds'),
                           greatest=c(TRUE))

# REDUCE
# drop all clusters with MAD scores less than .5
mad_scrs <- calc_mad(genus_only, genus_only@cids)
keep <- names(mad_scrs)[mad_scrs > .8]
genus_only <- drop_cls(genus_only, keep)

# SUMMARISE
smmry <- summary(genus_only)
smmry <- smmry[order(smmry$N_taxa, decreasing=TRUE), ]
write.csv(smmry, file.path('figures', 'best_clusters_primates.csv'))

# KEEP ONLY TOP TEN
slctd_smmry <- smmry[1:10, ]
slctd_smmry$ID <- as.numeric(slctd_smmry$ID)
slctd <- drop_cls(genus_only, as.character(slctd_smmry$ID))

# OUTPUT
# write out top 10 clusters with most taxa
sqfls <- list.files(wd, pattern='sequences[0-9]+.fasta')
rmFls(file.path(wd, sqfls))
for(i in seq_along(slctd@cids)) {
  cid <- slctd@cids[i]
  sids <- slctd@cls[[cid]]@sids
  txids <- get_txids(slctd, cid=cid, rnk='genus')
  scnms <- get_tx_slot(slctd, txids, 'scnm')
  n <- sapply(seq_along(scnms), function(x) 
    sum(scnms[x] == scnms[x:length(scnms)]))
  sq_nm <- paste0(scnms, '_', n)
  infile <- file.path(wd, paste0('sequences', i, '.fasta'))
  write_sqs(phylota=slctd, outfile=infile, sid=sids,
            sq_nm=sq_nm)
}

# ALIGN
alfls <- list.files(wd, pattern='alignment[0-9]+.fasta')
rmFls(file.path(wd, alfls))
for(i in seq_along(slctd@cids)) {
  inpt <- file.path(wd, paste0('sequences', i,
                               '.fasta'))
  otpt <- file.path(wd, paste0('alignment', i,'.fasta'))
  system(paste0('mafft --auto ', inpt, ' > ', otpt))
}
