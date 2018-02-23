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

getNNs <- function(phylo) {
  tree <- phylo2TM(phylo)
  res <- vector('list', length=tree['ntips'])
  names(res) <- tree@tips
  for(tp in tree@tips) {
    res[[tp]] <- treeman::getNdSstr(tree, id=tp)
  }
  res
}

calcNNs <- function(phylo) {
  # Proportion of tips who are sister with their genus
  nns <- getNNs(phylo)
  res <- sapply(seq_along(nns), function(x) {
    pttrn <- sub('_[0-9]+', '', names(nns)[[x]])
    any(grepl(pttrn, nns[[x]]))
  })
  sum(res)/length(res)
}

taxonomicallyInform <- function(tree, parent) {
  # Look up taxonomic ranks for nodes using GNR
  txnyms <- treeman::searchTxnyms(tree, cache=TRUE,
                                  infer=TRUE, clean=FALSE,
                                  parent=parent)
  tree <- treeman::setTxnyms(tree, txnyms)
  tree
}

reduceToFamily <- function(phylo, parent, tp_gls) {
  calc <- function(lng) {
    psslbs <- lng[lng %in% tp_gls]
    psslbs[length(psslbs)]
  }
  # convert phylo to TreeMan
  tree <- phylo2TM(phylo)
  tree <- taxonomicallyInform(tree, parent)
  lngs <- treeman::getNdsLng(tree, tree['tips'])
  gls <- unlist(sapply(lngs, calc))
  to_keep <- names(gls)[!duplicated(gls) & !is.na(gls)]
  nw_nms <- gls[!duplicated(gls) & !is.na(gls)]
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
