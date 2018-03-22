drpOvrhngs <- function(als, ctff=.75) {
  for(i in seq_along(als)) {
    gps_mtrx <- matrix(FALSE, ncol=nchar(als[[i]][[1]]),
                       nrow=length(als[[i]]))
    for(j in seq_along(als[[i]])) {
      pull <- gregexpr('-', als[[i]][[j]])[[1]]
      if(pull[1] == -1) {
        next
      }
      gps_mtrx[j, pull] <- TRUE
    }
    prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
    ovrlppng <- which(prp_mssng > ctff)
    strt <- ovrlppng[1]
    end <- ovrlppng[length(ovrlppng)]
    for(j in seq_along(als[[i]])) {
      als[[i]][[j]] <- substr(als[[i]][[j]], start=strt, stop=end)
    }
  }
  als
}

spltUp <- function(al, ctff=.75, mn_lngth=500) {
  gps_mtrx <- matrix(FALSE, ncol=nchar(al[[1]]),
                     nrow=length(al))
  for(i in seq_along(al)) {
    pull <- gregexpr('-', al[[i]])[[1]]
    gps_mtrx[i, pull] <- TRUE
  }
  prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
  strts <- which((prp_mssng[-1] > ctff) & (prp_mssng[-1*length(prp_mssng)] < ctff)) + 1
  ends <- which((prp_mssng[-1] < ctff) & (prp_mssng[-1*length(prp_mssng)] > ctff))
  lngths <- ends - strts
  pull <- lngths > mn_lngth
  strts <- strts[pull]
  ends <- ends[pull]
  new_als <- vector('list', length=length(strts))
  for(i in seq_along(strts)) {
    strt <- strts[i]
    end <- ends[i]
    new_al <- vector('list', length=length(al))
    for(j in seq_along(al)) {
      new_al[[j]] <- substr(al[[j]], start=strt, stop=end)
    }
    names(new_al) <- names(al)
    new_als[[i]] <- new_al
  }
  new_als
}


readSqs <- function(fl) {
  all_data <- readLines(fl)
  sqs <- list()
  for(i in seq_along(all_data)) {
    bit <- all_data[[i]]
    if(grepl(pattern='^>', x=bit)) {
      nm <- sub(pattern='^>', '', x=bit)
      sqs[[nm]] <- ''
    } else {
      sqs[[nm]] <- paste0(sqs[[nm]], bit)
    }
  }
  sqs
}

writeSqs <- function(sqs, fl) {
  fasta <- ''
  for(i in seq_along(sqs)) {
    sq <- sqs[[i]]
    dfln <- names(sqs)[[i]]
    fasta <- paste0(fasta, '>', dfln, '\n',
                    sq, '\n\n')
  }
  cat(fasta, file=fl)
}

partition <- function(lngths, fl) {
  gene <- strt <- 1
  prttn_txt <- ''
  for(lngth in lngths) {
    end <- lngth + strt - 1
    prttn_txt <- paste0(prttn_txt, 'DNA, gene',
                        gene, ' = ', strt, '-',
                        end, '\n')
    strt <- end + 1
    gene <- gene + 1
  }
  cat(prttn_txt, file=fl)
}