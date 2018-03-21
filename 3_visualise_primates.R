# Visualise the primates tree

# LIBS
library(ape)
library(doMC)
source(file.path('tools', 'treeman_tools.R'))

# VARS
wd <- 'primates'
prosimians <- read.delim(file.path('expected', 'prosimians_list.txt'),
                         stringsAsFactors=FALSE)[ ,1]

# INPUT
trstr <- readLines(file.path(wd, 'consensus.tre'))
trstr <- gsub(':[0-9.]+', '', trstr)
trstr <- gsub('(\\[|\\])', '', trstr)
tree <- read.tree(text=trstr)
expctd <- read.nexus(file.path('expected', 'primates.nex'))

# REDUCE TREE
nns_dst <- calcNNs(tree)
tp_lbls <- tree[['tip.label']]
tp_lbls <- sub('_[0-9]+', '', tp_lbls)
to_drp <- tree[['tip.label']][duplicated(tp_lbls)]
tree <- drop.tip(tree, to_drp)
tree[['tip.label']] <- sub('_[0-9]+', '', tree[['tip.label']])
tree[['tip.label']] <- sub('_.*', '', tree[['tip.label']])

# REDUCE EXPECTED TO GENUS-LEVEL
genus_labels <- sub('_.*$', '', expctd[['tip.label']])
to_drp <- expctd$tip.label[duplicated(genus_labels)]
expctd <- drop.tip(expctd, tip=to_drp)
expctd$tip.label <- sub('_.*$', '', expctd$tip.label)

# REROOT
prosimians <- prosimians[prosimians %in% tree$tip.label]
tree <- unroot(tree)
tree <- root(tree, outgroup=prosimians, resolve.root=TRUE)

# HOW MANY GENERA?
esrch <- rentrez::entrez_search(db='taxonomy', term='txid9443[Subtree] AND genus[Rank]')
genera_counts <- data.frame('actual'=esrch[['count']],
                            'obs'=length(tree$tip.label))
genera_counts[['prp']] <- genera_counts[['obs']]/genera_counts[['actual']]
write.csv(genera_counts, file=file.path('results', 'primates_counts.csv'),
          row.names=FALSE)

# BRANCH SUPPORTS
spprt <- tree$node.label
spprt <- suppressWarnings(as.numeric(spprt))
spprt[is.na(spprt)] <- 0
nd_lbls <- rep('', length(spprt))
nd_lbls[spprt > 50] <- '*'
nd_lbls[spprt > 75] <- '**'
nd_lbls[spprt > 95] <- '***'

# PLOT
png(file.path('results', 'primates.png'), width=2000, height=2000)
par(mar=c(.1, .1, .1, .1))
plot(tree, edge.width=4, cex=2)
nodelabels(text=nd_lbls, frame='none', cex=2.5, adj=-.25)
dev.off()
saveRDS(tree, file.path('figures', 'primates_tree.RData'))

# DROP UNSHARED TIPS
to_drp <- expctd$tip.label[!expctd$tip.label %in% tree$tip.label]
expctd <- drop.tip(expctd, tip=to_drp)
to_drp <- tree$tip.label[!tree$tip.label %in% expctd$tip.label]
tree_cmp <- drop.tip(tree, tip=to_drp)
tree_cmp <- ladderize(tree_cmp, right=TRUE)
expctd <- ladderize(expctd, right=TRUE)

# DISTS
dsts <- compareTrees(tree_cmp, expctd)
write.csv(dsts, file=file.path('results', 'primates_dst.csv'),
          row.names=FALSE)

# REDUCE
tree_cmp <- reduceToFamily(tree_cmp, parent='Primates', tp_gls='family')
expected <- reduceToFamily(expected, parent='Primates')

# COPLOT
if(!exists('dsts')) {
  dsts <- read.csv(file.path('results', 'primates_dst.csv'))
}
if(!exists('genera_counts')) {
  genera_counts <- read.csv(file.path('results', 'primates_counts.csv'))
}
png(file.path('results', 'primates_coplot.png'), width=2000, height=2000)
par(cex=2, mar=c(1,.1,.1,.1))
suppressWarnings(cophyloplot(tree_cmp, expctd, space=10,
                             gap=5))
mtext(text='phylotaR', side=3, cex=2.5, line=-1.5, adj=.25)
mtext(text='Expected', side=3, cex=2.5, line=-1.5, adj=.75)
mtext(text=paste0('RF ', round(dsts[['rf_dst']], 3),
                  ' | TRP ', round(dsts[['trp_dst']], 3)),
      side=1, line=-1, adj=.1, cex=2.5)
mtext(text=paste0('Obs N. ', genera_counts[['obs']],
                  ' | Actual N. ', genera_counts[['actual']],
                  ' (', signif(genera_counts[['prp']], 2), ')'),
      side=1, line=-1, adj=.9, cex=2.5)
dev.off()
