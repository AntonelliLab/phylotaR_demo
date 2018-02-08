# TREE MANIPULATION FUNCTIONS
# Requires different package with different tree structure

# FUNCTIONS
tm2Phylo <- function(tree) {
  treeman::writeTree(tree, file='temp.tre')
  phylo <- read.tree(file='temp.tre')
  file.remove('temp.tre')
  phylo
}

phylo2TM <- function(phylo) {
  write.tree(phylo, file='temp.tre')
  tree <- treeman::readTree(file='temp.tre')
  file.remove('temp.tre')
  tree
}

taxonomicallyInform <- function(tree, parent) {
  # Look up taxonomic ranks for nodes using GNR
  txnyms <- treeman::searchTxnyms(tree, cache=TRUE, infer=TRUE,
                         clean=TRUE, parent=parent)
  tree <- treeman::setTxnyms(tree, txnyms)
  tree
}

reduceToFamily <- function(phylo, parent) {
  # convert phylo to TreeMan
  tree <- phylo2TM(phylo)
  tree <- taxonomicallyInform(tree, parent)
  lngs <- treeman::getNdsLng(tree, tree['tips'])
  # can't use named ranks, choosing level up
  fmls <- vapply(lngs, function(x) x[[length(x)-1]], '')
  # ignore all ambiguous entries
  fmls <- fmls[fmls != '']
  to_keep <- names(fmls)[!duplicated(fmls)]
  nw_nms <- fmls[!duplicated(fmls)]
  to_drop <- tree['tips'][!tree['tips'] %in% to_keep]
  tree <- treeman::rmTips(tree, to_drop)
  tree <- treeman::setNdsID(tree, to_keep, nw_nms)
  # return as phylo
  tm2Phylo(tree)
}

compareTrees <- function(phylo_1, phylo_2, parallel=FALSE) {
  tree_1 <- phylo2TM(phylo_1)
  tree_2 <- phylo2TM(phylo_2)
  tree_1 <- suppressMessages(treeman::addNdmtrx(tree_1))
  tree_2 <- suppressMessages(treeman::addNdmtrx(tree_2))
  rf_dst <- treeman::calcDstRF(tree_1, tree_2)
  trp_dst <- treeman::calcDstTrp(tree_1, tree_2, parallel=parallel)
  data.frame('rf_dst'=rf_dst, 'trp_dst'=trp_dst)
}
