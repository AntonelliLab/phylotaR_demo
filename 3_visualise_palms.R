# Visualise the palms tree

# LIBS
library(ape)
source(file.path('tools', 'treeman_tools.R'))

# VARS
wd <- 'palms'
calamoideae <- read.delim(file.path('expected', 'calamoideae_list.txt'),
                          stringsAsFactors=FALSE)[ ,1]

# INPUT
trstr <- readLines(file.path(wd, 'consensus.tre'))
trstr <- gsub(':[0-9.]+', '', trstr)
trstr <- gsub('(\\[|\\])', '', trstr)
tree <- read.tree(text=trstr)
expctd <- read.nexus(file.path('expected', 'palms.nex'))
expctd[['tip.label']] <- sub('_.*', '', expctd[['tip.label']])

# REDUCE TREE
nns_dst <- calcNNs(tree)
tp_lbls <- tree[['tip.label']]
tp_lbls <- sub('_[0-9]+', '', tp_lbls)
to_drp <- tree[['tip.label']][duplicated(tp_lbls)]
tree <- drop.tip(tree, to_drp)
tree[['tip.label']] <- sub('_[0-9]+', '', tree[['tip.label']])
tree[['tip.label']] <- sub('_.*', '', tree[['tip.label']])

# REDUCE EXPECTED TO TRIBE-LEVEL
tp_gls <- tree[['tip.label']]
expctd <- reduceToFamily(expctd, parent='Arecaceae',
                         tp_gls=tp_gls)

# REROOT
calamoideae <- calamoideae[calamoideae %in% tree$tip.label]
tree <- unroot(tree)
tree <- root(tree, outgroup=calamoideae, resolve.root=TRUE)
write.tree(tree, file=file.path('results', 'palms.tre'))

# HOW MANY TRIBES?
esrch <- rentrez::entrez_search(db='taxonomy', term='txid4710[Subtree] AND tribe[Rank]')
genera_counts <- data.frame('actual'=esrch[['count']],
                            'obs'=length(tree$tip.label))
genera_counts[['prp']] <- genera_counts[['obs']]/genera_counts[['actual']]
write.csv(genera_counts, file=file.path('results', 'palms_counts.csv'),
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
png(file.path('results', 'palms.png'), width=2000, height=2000)
par(mar=c(.1, .1, .1, .1))
plot(tree, edge.width=4, cex=3.5)
nodelabels(text=nd_lbls, frame='none', cex=2.5, adj=-.25)
dev.off()
tree$node.label <- nd_lbls
saveRDS(tree, file.path('figures', 'palms_tree.RData'))

# DROP UNSHARED TIPS
to_drp <- tree$tip.label[!tree$tip.label %in% expctd$tip.label]
tree_cmp <- drop.tip(tree, tip=to_drp)

tree_cmp <- ladderize(tree_cmp)
expctd <- ladderize(expctd)

# DISTS
dsts <- compareTrees(tree_cmp, expctd)
write.csv(dsts, file=file.path('results', 'palms_dst.csv'),
          row.names=FALSE)

# COPLOT
if(!exists('dsts')) {
  dsts <- read.csv(file.path('results', 'palms_dst.csv'))
}
if(!exists('genera_counts')) {
  genera_counts <- read.csv(file.path('results', 'palms_counts.csv'))
}

assoc <- cbind(tree_cmp$tip.label, tree_cmp$tip.label)
tree_cmp <- rotate(tree_cmp, c("Caryoteae", "Borasseae"))
tree_cmp <- rotate(tree_cmp, c("Calameae", "Lepidocaryeae"))
tree_cmp <- rotate(tree_cmp, c("Euterpeae", "Pelagodoxeae"))

png(file.path('results', 'palms_coplot.png'), width=2000, height=2000)
par(cex=2, mar=c(1,.1,.1,.1))
suppressWarnings(cophyloplot(tree_cmp, expctd, space=80,
                             gap=5, edge.width=2, assoc=assoc))
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
