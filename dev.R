
# rnk <- 'tribe'
# keep_hghr_taxa <- TRUE
# which.max(n_taxa)
# id <- clstrs_obj@clstr_ids[[1]]
# mn_pambg <- 0.5

fltrClstrSqs <- function(clstrs_obj, id, rnk='species',
                         n_ech=1, mn_pambg=0.1, hghr_tx=TRUE) {
  slct <- function(txid) {
    pssbls <- sqlns[which(txid == txids)]
    mx_n <- ifelse(length(pssbls) > n_ech, n_ech, length(pssbls))
    pssbls <- sort(x=pssbls, decreasing=TRUE)[1:mx_n]
    names(pssbls)
  }
  clstr <- clstrs_obj@clstrs[[id]]
  pambgs <- getSqAmbs(clstrs_obj=clstrs_obj, id=id)
  gis <- names(pambgs)[pambgs < mn_pambg]
  txids <- getSqTx(clstrs_obj, sid=gis, rnk=rnk,
                   hghr_tx=hghr_tx)
  sqlns <- getSqLns(clstrs_obj, sid=gis)
  names(sqlns) <- gis
  unqids <- unique(txids)
  keep <- unlist(lapply(unqids, slct))
  drpSqs(clstrs_obj=clstrs_obj, cid=id, sid=keep)
}


getSqTx <- function(clstrs_obj, cid=NULL, sid=NULL,
                    rnk=FALSE, hghr_tx=TRUE) {
  if(!is.null(cid)) {
    txids <- clstrs_obj@clstrs[[cid]][['tis']]
  } else {
    txids <- sapply(sid,
                    function(x) clstrs_obj@sqs[[x]][['ti']])
  }
  if(rnk != FALSE) {
    txids <- getIDFrmTxdct(txdct=clstrs_obj@txdct,
                           id=txids, rnk=rnk,
                           hghr_tx=hghr_tx)
  }
  txids
}

getIDFrmTxdct <- function(txdct, id, ret=c('TaxId', 'ScientificName'),
                          rnk=c("superkingdom", "kingdom", "phylum",
                                 "subphylum", "class", "superorder",
                                 "order", "suborder", "infraorder",
                                 "parvorder", "family", "subfamily",
                                 "tribe", "subtribe", "genus",
                                 "species", "subspecies"),
                          hghr_tx=TRUE) {
  calc <- function(id) {
    rcrd <- txdct[[id]]
    rnks <- sapply(rcrd[['LineageEx']],
                   function(x) x[['Rank']])
    ids <- sapply(rcrd[['LineageEx']],
                  function(x) x[[ret]])
    mtch <- which(rnk == rnks)
    if(length(mtch) == 0) {
      # if no match, use lowest available rnk
      mtch <- length(rnks)
    }
    ids[[mtch[[1]]]]
  }
  rnk <- match.arg(rnk)
  ret <- match.arg(ret)
  sapply(as.character(id), calc)
}

rmSqs <- function(clstrs_obj, blcklst) {
  for(cid in clstrs_obj@clstr_ids) {
    clstr <- clstrs_obj@clstrs[[cid]]
    if(any(as.character(blcklst) %in% as.character(clstr$gis))) {
      pull <- !(as.character(clstr$gis) %in% as.character(blcklst))
      clstr[['tis']] <- clstr[['tis']][pull]
      clstr[['gis']] <- clstr[['gis']][pull]
      clstr[['n_gi']] <- clstr$n_gi - 1
      clstr[['n_ti']] <- length(unique(clstr[['tis']]))
      clstrs_obj@clstrs[[cid]] <- clstr
    }
  }
  clstrs_obj
}
