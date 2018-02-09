# Root and visualise the palms tree

# LIBS
library(ape)
library(doMC)
source(file.path('tools', 'treeman_tools.R'))

# VARS
wd <- 'palms'
calamoideae <- c('Calamus', 'Ceratolobus', 'Daemonorops',
                 'Pogonotium', 'Retispatha', 'Korthalsia',
                 'Metroxylon', 'Pigafetta', 'Myrialepis',
                 'Plectocomia', 'Plectocomiopsis', 'Eleiodoxa',
                 'Salacca', 'Eugeissona', 'Eremospatha',
                 'Laccosperma', 'Oncocalamus', 'Lepidocaryum',
                 'Mauritia', 'Mauritiella', 'Raphia')
# calculating trip dist takes time, run in parallel
# use doSnow for Windows
registerDoMC(cores=4)

# INPUT
tree <- read.tree(file=file.path(wd, 'tree.tre'))
expctd <- read.nexus(file.path('expected', 'palms.nex'))

# REROOT
calamoideae <- calamoideae[calamoideae %in% tree$tip.label]
tree <- unroot(tree)
tree <- root(tree, outgroup=calamoideae, resolve.root=TRUE)

# HOW MANY GENERA?
esrch <- rentrez::entrez_search(db='taxonomy', term='txid4710[Subtree] AND genus[Rank]')
genera_counts <- data.frame('actual'=esrch[['count']],
                            'obs'=length(tree$tip.label))
genera_counts[['prp']] <- genera_counts[['obs']]/genera_counts[['actual']]
write.csv(genera_counts, file=file.path('figures', 'palms_counts.csv'),
          row.names=FALSE)

# PLOT
png(file.path('figures', 'palms.png'), width=2000, height=2000)
par(mar=c(.1, .1, .1, .1))
plot(tree, type='radial', cex=1.5)
dev.off()

# DROP UNSHARED TIPS
expctd$tip.label <- sub('_.*$', '', expctd$tip.label)
to_drp <- expctd$tip.label[!expctd$tip.label %in% tree$tip.label]
expctd <- drop.tip(expctd, tip=to_drp)
to_drp <- tree$tip.label[!tree$tip.label %in% expctd$tip.label]
tree <- drop.tip(tree, tip=to_drp)

# REDUCE TO FAMILY-LEVEL
tree_fmly <- reduceToFamily(tree, parent='Arecaceae')
expctd_fmly <- reduceToFamily(expctd, parent='Arecaceae')

# DISTS
dsts <- compareTrees(tree, expctd, parallel=TRUE)
write.csv(dsts, file=file.path('figures', 'palms_dst.csv'),
          row.names=FALSE)

# COPLOT
if(!exists('dsts')) {
  dsts <- read.csv(file.path('figures', 'palms_dst.csv'))
}
if(!exists('genera_counts')) {
  genera_counts <- read.csv(file.path('figures', 'palms_counts.csv'))
}
pull <- tree_fmly$tip.label %in% expctd_fmly$tip.label
assoc <- cbind(tree_fmly$tip.label[pull], tree_fmly$tip.label[pull])
png(file.path('figures', 'palms_coplot.png'), width=2000, height=2000)
par(cex=2, mar=c(1,.1,.1,.1))
suppressWarnings(cophyloplot(tree_fmly, expctd_fmly, space=10,
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
