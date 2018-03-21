# Identify best clusters

# LIBS
library(ggplot2)
library(phylotaR)
source(file.path('tools', 'selection_tools.R'))

# VARS
wd <- 'palms'

# INPUT
all_cls <- read_phylota(wd)

# CLTYPE STATS
# drop clusters of 10
n_taxa <- get_cl_slot(all_cls, cid=all_cls@cids, slt_nm='ntx')
keep <- names(n_taxa)[n_taxa > 10]
all_cls <- drop_cls(all_cls, keep)
table(sapply(all_cls@cls@cls, function(x) x@typ))

# PLOT
# tribe treemap
tribe_txids <- unique(get_txids(phylota=all_cls,
                                txids=all_cls@txids,
                                rnk='tribe'))
tribe_txids <- tribe_txids[tribe_txids != '']
txnms <- get_tx_slot(phylota=all_cls, txid=tribe_txids,
                     slt_nm='scnm')
txnms <- sort(txnms, decreasing=TRUE)
p <- plot_phylota_treemap(phylota=all_cls, txids=tribe_txids,
                          txnms=txnms, area='nsq', fill='ncl')
png(file.path('results', 'palms_tx_treemap.png'), width=2000, height=2000)
print(p + theme(legend.position='none'))
dev.off()
saveRDS(p, file.path('figures', 'palms_tx_treemap.RData'))
# cluster treemap
p <- plot_phylota_treemap(phylota=all_cls, cids=all_cls@cids,
                          area='nsq', fill='ntx')
png(file.path('results', 'palms_cl_treemap.png'), width=2000, height=2000)
print(p + theme(legend.position='none'))
dev.off()
saveRDS(p, file.path('figures', 'palms_cl_treemap.RData'))

# REDUCE TO TRIBE
tribe_only <- drop_by_rank(all_cls, rnk='tribe', n=2,
                           choose_by= c("pambgs", "age",
                                        "nncltds"),
                           greatest = c(FALSE, FALSE,
                                        TRUE))

# SUMMARISE AND FILTER
# count n genera per cid
n_tribes <- sapply(tribe_only@cids, function(x) {
  length(unique(get_txids(phylota=tribe_only, cid=x, rnk='tribe')))
})
smmry <- summary(tribe_only)
smmry[['N_taxa']] <- n_tribes
# drop all with fewer than 0.5 MAD
smmry <- smmry[smmry[['MAD']] > 0.5, ]
smmry <- smmry[order(smmry$N_taxa, decreasing=TRUE), ]

# SELECT
slctd_smmry <- smmry[1:10, ]
slctd_smmry$ID <- as.numeric(slctd_smmry$ID)
slctd <- drop_cls(tribe_only, as.character(slctd_smmry$ID))
write.csv(slctd_smmry, file.path('results', 'best_clusters_palms.csv'))

# PLOT
tribes <- unique(get_txids(phylota=all_cls,
                           txids=all_cls@txids, rnk='tribe'))
tribes <- tribes[tribes != '']
txnms <- get_tx_slot(phylota=slctd, txid=tribes, slt_nm='scnm')
ordd <- order(txnms, decreasing=TRUE)
p <- plot_phylota_pa(phylota=slctd, cids=slctd@cids,
                     txids=tribes[ordd], txnms=txnms[ordd])
png(file.path('results', 'palms_cl_pamap.png'), width=2000, height=2000)
print(p)
dev.off()
saveRDS(p, file.path('figures', 'palms_cl_pamap.RData'))

# OUTPUT
# write out top 10 clusters with most taxa
sqfls <- list.files(wd, pattern='sequences[0-9]+.fasta')
rmFls(file.path(wd, sqfls))
for(i in seq_along(slctd@cids)) {
  cid <- slctd@cids[i]
  sids <- slctd@cls[[cid]]@sids
  txids <- get_txids(slctd, cid=cid, rnk='tribe')
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
