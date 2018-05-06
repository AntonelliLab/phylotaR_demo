# Identify best clusters

# LIBS
library(ggplot2)
library(phylotaR)
source(file.path('tools', 'selection_tools.R'))

# VARS
wd <- 'primates'

# INPUT
all_cls <- read_phylota(wd)

# RECREATE PhyLoTa TABLE
# get phylogenetically inform. clusters
n_taxa <- get_cl_slot(all_cls, cid = all_cls@cids, slt_nm = 'ntx')
cltyps <- get_cl_slot(all_cls, all_cls@cids, slt_nm = 'typ')
#picls <- cltyps %in% c('subtree', 'direct') & n_taxa >= 4
picls <- n_taxa >= 4
picls <- all_cls@cids[picls]
# build table
prnts <- get_cl_slot(all_cls, cid = picls, slt_nm = 'prnt')
n_taxa <- get_cl_slot(all_cls, cid = picls, slt_nm = 'ntx')
n_sqs <- get_cl_slot(all_cls, cid = picls, slt_nm = 'nsqs')
n_gnra <- sapply(picls, function(x) {
  length(unique(get_txids(all_cls, cid = x, rnk = 'genus')))
})
mxsqlngs <- sapply(picls, function(x) {
  max(get_sq_slot(all_cls, cid = x, slt_nm = 'nncltds'))
})
mnsqlngs <- sapply(picls, function(x) {
  min(get_sq_slot(all_cls, cid = x, slt_nm = 'nncltds'))
})
lngst <- sapply(picls, function(x) {
  sids <- all_cls[[x]]@sids
  lngths <- get_sq_slot(all_cls, sid = sids, slt_nm = 'nncltds')
  sids[which.max(lngths)]
})
deflns <- get_sq_slot(all_cls, sid = lngst, slt_nm = 'dfln')
mads <- sapply(picls, calc_mad, phylota = all_cls)
phylota_table <- data.frame('Cluster.ID' = picls, 'Parent' = prnts,
                            'TaxIDs' = n_taxa, 'GIs' = n_sqs,
                            'Genera' = n_gnra, 'L.min' = mnsqlngs,
                            'L.max' = mxsqlngs, 'MAD' = mads,
                            'Defline.of.longest.sequence' = deflns)
ordr <- order(phylota_table[['TaxIDs']], decreasing = TRUE)
phylota_table <- phylota_table[ordr, ]
write.csv(x = phylota_table, file = file.path('results', 'primates_phylota.csv'),
          row.names = FALSE)

# CLTYPE STATS
# drop clusters of 10
n_taxa <- get_cl_slot(all_cls, cid=all_cls@cids, slt_nm='ntx')
keep <- names(n_taxa)[n_taxa > 10]
all_cls <- drop_cls(all_cls, keep)
table(sapply(all_cls@cls@cls, function(x) x@typ))

# PLOT
# genus treemap
genus_txids <- unique(get_txids(phylota=all_cls,
                                txids=all_cls@txids,
                                rnk='genus'))
genus_txids <- genus_txids[genus_txids != '']
txnms <- get_tx_slot(phylota=all_cls, txid=genus_txids,
                     slt_nm='scnm')
txnms <- sort(txnms, decreasing=TRUE)
p <- plot_phylota_treemap(phylota=all_cls, txids=genus_txids,
                          txnms=txnms, area='nsq', fill='ncl')
png(file.path('results', 'primates_tx_treemap.png'), width=2000, height=2000)
print(p + theme(legend.position='none'))
dev.off()
saveRDS(p, file.path('figures', 'primates_tx_treemap.RData'))
# cluster treemap
p <- plot_phylota_treemap(phylota=all_cls, cids=all_cls@cids,
                          area='nsq', fill='ntx')
png(file.path('results', 'primates_cl_treemap.png'), width=2000, height=2000)
print(p + theme(legend.position='none'))
dev.off()
saveRDS(p, file.path('figures', 'primates_cl_treemap.RData'))

# REDUCE TO GENUS
genus_only <- drop_by_rank(all_cls, rnk='genus', n=2,
                           choose_by= c("nncltds", 'pambgs'),
                           greatest = c(TRUE, FALSE))

# SUMMARISE AND FILTER
# count n genera per cid
n_genera <- sapply(genus_only@cids, function(x) {
  length(unique(get_txids(phylota=genus_only, cid=x, rnk='genus')))
})
smmry <- summary(genus_only)
smmry[['N_taxa']] <- n_genera
smmry <- smmry[smmry[['MAD']] > 0.75, ]
smmry <- smmry[order(smmry$N_taxa, decreasing=TRUE), ]

# SELECT
slctd_smmry <- smmry[1:10, ]
slctd_smmry$ID <- as.numeric(slctd_smmry$ID)
slctd <- drop_cls(genus_only, as.character(slctd_smmry$ID))
write.csv(slctd_smmry, file.path('results', 'best_clusters_primates.csv'))

# PLOT
genera <- unique(get_txids(phylota=all_cls,
                           txids=all_cls@txids, rnk='genus'))
genera <- genera[genera != '']
txnms <- get_tx_slot(phylota=slctd, txid=genera, slt_nm='scnm')
ordd <- order(txnms, decreasing=TRUE)
p <- plot_phylota_pa(phylota=slctd, cids=slctd@cids,
                     txids=genera[ordd], txnms=txnms[ordd])
png(file.path('results', 'primates_cl_pamap.png'), width=2000, height=2000)
print(p)
dev.off()
saveRDS(p, file.path('figures', 'primates_cl_pamap.RData'))

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
