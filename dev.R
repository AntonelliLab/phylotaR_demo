mkdata <- function(cid) {
  cl <- phylota@cls@cls[[which(cid == phylota@cids)]]
  value <- apply(X=tis_mtrx[ ,cl@sids], 1, any)
  value <- as.numeric(value)
  data.frame(txid=as.character(txids),
             cid=cid, value=value)
}
value <- NULL
# gen p_data
tis_mtrx <- phylotaR:::mk_txid_in_sq_mtrx(phylota=phylota,
                               txids=txids)
p_data <- lapply(cids, mkdata)
p_data <- do.call('rbind', p_data)
p_data[['cnm']] <- 
  cnms[match(p_data[['cid']], cids)]
p_data[['txnm']] <- 
  txnms[match(p_data[['txid']], txids)]
# reorder
p_data[['cnm']] <- factor(p_data[['cnm']],
                          levels=cnms,
                          ordered=TRUE)
p_data[['txnm']] <- factor(p_data[['txnm']],
                           levels=txnms,
                           ordered=TRUE)
# plot
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
cnm <- txnm <- NULL
ggplot2::ggplot(p_data, ggplot2::aes(cnm, txnm)) +
  ggplot2::geom_tile(ggplot2::aes(fill=value)) +
  ggplot2::xlab('') + ggplot2::ylab('') +
  ggplot2::scale_fill_gradient(low='#e0e0e0',
                               high='#303030') +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position='none')
