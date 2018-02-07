library(phylotaR)
wd <- '/home/dom/Desktop/testing_phylotaR/primates'

# IDENTIFY BEST CLUSTER
# generate cluster object
clstrs_obj <- genClstrsObj(wd)
# get n taxa per cluster
n_taxa <- getClNTx(clstrs_obj,
                   clstrs_obj@clstr_ids)
# drop all clusters with fewer than 50 taxa
keep <- names(n_taxa)[n_taxa > 100]
clstrs_obj <- drpClstrs(clstrs_obj, keep)
# filter each cluster in the cluster object
#  - keep only sequences repres. genus or higher
#  - drop all sequences with % ambiguous bases of more than 50%
#  - keep only the largest sequence of each genus
for(id in clstrs_obj@clstr_ids) {
  clstrs_obj<- fltrClstrSqs(clstrs_obj, id=id,
                            rank='genus', mn_pambg=0.5)
}
# drop all clusters with MAD scores less than .95
mad_scrs <- getClMAD(clstrs_obj, clstrs_obj@clstr_ids)
keep <- names(mad_scrs)[mad_scrs > .95]
clstrs_obj <- drpClstrs(clstrs_obj, keep)
# look at the median sequence length of each cluster
(sqlngs <- getClMdLn(clstrs_obj, clstrs_obj@clstr_ids))
# check n taxa
(n_taxa <- getClNTx(clstrs_obj, clstrs_obj@clstr_ids))
# drop cluster with less than max sequences
keep <- names(n_taxa)[n_taxa == max(n_taxa)]
clstrs_obj <- drpClstrs(clstrs_obj, keep)
# look at most common words in definition lines of each cluster
for(id in clstrs_obj@clstr_ids) {
  print(getSqDfs(clstrs_obj, id=id, prse=0.2))
}
# select the one with the longest seq
(sqlngs <- getClMdLn(clstrs_obj, clstrs_obj@clstr_ids))
# select a cluster ID
cid <- clstrs_obj@clstr_ids[[which.max(sqlngs)]]
getSqDfs(clstrs_obj, id=cid, prse=0.1)
# get its sequences IDs
sids <- clstrs_obj@clstrs[[cid]][['gis']]
# get genus names for fasta def. lines
scnms <- getIDFrmTxdct(clstrs_obj@txdct, ret='ScientificName',
                       id=names(sids), rank='genus')
infile <- file.path(wd, 'sequences.fasta')
outfile <- file.path(wd, 'alignment.fasta')
writeSqs(sids=sids, dflns=scnms, flpth=infile)

# mafft --auto sequences.fasta > alignment.fasta
# raxmlHPC -m GTRGAMMA -f a -N 10 -p 1234 -x 1234 -n primates -s alignment.fasta