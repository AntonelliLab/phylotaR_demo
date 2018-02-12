# Root and visualise the palms tree
# ... also calc tree dists

# LIBS
library(ape)
source(file.path('tools', 'treeman_tools.R'))

# VARS
wd <- 'primates'
prosimians <- read.delim(file.path(wd, 'prosimians'),
                         stringsAsFactors=FALSE)[ ,1]
prosimians <- c('Hapalemur', 'Lemur', 'Varecia', 'Eulemur',
                'Avahi', 'Propithecus', 'Mirza', 'Daubentonia',
                'Galago', 'Arctocebus', 'Perodicticus', 'Loris',
                'Nycticebus', 'Galagoides', 'Euoticus')

# INPUT
trstr <- readLines(file.path(wd, 'consensus.tre'))
trstr <- gsub(':[0-9.]+', '', trstr)
trstr <- gsub('(\\[|\\])', '', trstr)
tree <- read.tree(text=trstr)
expctd <- read.nexus(file.path('expected', 'primates.nex'))

# REROOT
prosimians <- prosimians[prosimians %in% tree$tip.label]
tree <- unroot(tree)
tree <- root(tree, outgroup=prosimians, resolve.root=TRUE)

# HOW MANY GENERA?
esrch <- rentrez::entrez_search(db='taxonomy', term='txid9443[Subtree] AND genus[Rank]')
genera_counts <- data.frame('actual'=esrch[['count']],
                            'obs'=length(tree$tip.label))
genera_counts[['prp']] <- genera_counts[['obs']]/genera_counts[['actual']]
write.csv(genera_counts, file=file.path('figures', 'primates_counts.csv'),
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
png(file.path('figures', 'primates.png'), width=2000, height=2000)
par(mar=c(.1, .1, .1, .1))
plot(tree, cex=1.5)
nodelabels(text=nd_lbls, frame='none', cex=1.5, adj=-1)
dev.off()

# REDUCE TO EXPCTD TO GENUS
genus_labels <- sub('_.*$', '', expctd$tip.label)
to_drp <- expctd$tip.label[duplicated(genus_labels)]
expctd <- drop.tip(expctd, tip=to_drp)
expctd$tip.label <- sub('_.*$', '', expctd$tip.label)

# DROP UNSHARED TIPS
to_drp <- expctd$tip.label[!expctd$tip.label %in% tree$tip.label]
expctd <- drop.tip(expctd, tip=to_drp)
to_drp <- tree$tip.label[!tree$tip.label %in% expctd$tip.label]
tree <- drop.tip(tree, tip=to_drp)

# DISTS
dsts <- compareTrees(tree, expctd)
write.csv(dsts, file=file.path('figures', 'primates_dst.csv'),
          row.names=FALSE)

# COPLOT
if(!exists('dsts')) {
  dsts <- read.csv(file.path('figures', 'primates_dst.csv'))
}
if(!exists('genera_counts')) {
  genera_counts <- read.csv(file.path('figures', 'primates_counts.csv'))
}
pull <- tree$tip.label %in% expctd$tip.label
assoc <- cbind(tree$tip.label[pull], tree$tip.label[pull])
png(file.path('figures', 'primates_coplot.png'), width=2000, height=2000)
par(cex=2, mar=c(1,.1,.1,.1))
suppressWarnings(cophyloplot(tree, expctd, space=10,
                             gap=5))
mtext(text='phylotaR', side=3, cex=2.5, line=-1.5, adj=.25)
mtext(text='Expected', side=3, cex=2.5, line=-1.5, adj=.75)
mtext(text=paste0('RF dist: ', dsts[['rf_dst']],
                  ' | TRP dist: ', dsts[['trp_dst']]),
      side=1, line=-1, adj=.1, cex=2.5)
mtext(text=paste0('Obs N. ', genera_counts[['obs']],
                  ' | Actual N. ', genera_counts[['actual']],
                  ' (', signif(genera_counts[['prp']], 2), ')'),
      side=1, line=-1, adj=.9, cex=2.5)
dev.off()