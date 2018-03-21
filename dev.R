  mkdata <- function(cid) {
    cl <- phylota@cls@cls[[which(cid == phylota@cids)]]
    value <- apply(X=txids_sqs[ ,cl@sids], 1, any)
    value <- factor(as.numeric(value))
    data.frame(txid=as.character(txids),
               cid=cid, value=value)
  }
  txidinsqs <- function(txid) {
    vapply(phylota@sids, function(x) {
      is_txid_in_sq(phylota=phylota, txid=txid, sid=x)
      }, logical(1))
  }
  txids_sqs <- lapply(txids, txidinsqs)
  txids_sqs <- do.call('rbind', txids_sqs)
  value <- NULL
  # gen p_data
  p_data <- lapply(cids, mkdata)
  p_data <- do.call('rbind', p_data)
  p_data[['cnm']] <- 
    cnms[match(p_data[['cid']], cids)]
  p_data[['txnm']] <- 
    txnms[match(p_data[['txid']], txids)]
  # reorder
  p_data[['cnm']] <- factor(p_data[['cnm']],
                            levels=cord,
                            ordered=TRUE)
  p_data[['txnm']] <- factor(p_data[['txnm']],
                             levels=txord,
                             ordered=TRUE)
  cnm <- txnm <- NULL
  ggplot2::ggplot(p_data, ggplot2::aes(cnm, txnm)) +
    ggplot2::geom_tile(ggplot2::aes(fill=value)) +
    ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::scale_fill_manual(values=c('white', 'black')) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position='none')
