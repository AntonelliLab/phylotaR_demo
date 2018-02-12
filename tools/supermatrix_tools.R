drpOvrhngs <- function(sqs, ctff=.75) {
  for(sq in sqs) {
    gps_mtrx <- matrix(FALSE, ncol=nchar(sq[[1]]),
                       nrow=length(sq))
    for(i in seq_along(sq)) {
      pull <- gregexpr('-', sq[[i]])[[1]]
      gps_mtrx[i, pull] <- TRUE
    }
    prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
    ovrlppng <- which(prp_mssng > ctff)
    strt <- ovrlppng[1]
    end <- ovrlppng[length(ovrlppng)]
    for(i in seq_along(sq)) {
      sq[[i]] <- substr(sq[[i]], start=strt, stop=end)
    }
  }
  sqs
}

spltUp <- function(sq, ctff=.75, mn_lngth=500) {
  gps_mtrx <- matrix(FALSE, ncol=nchar(sq[[1]]),
                     nrow=length(sq))
  for(i in seq_along(sq)) {
    pull <- gregexpr('-', sq[[i]])[[1]]
    gps_mtrx[i, pull] <- TRUE
  }
  prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
  strts <- which((prp_mssng[-1] > ctff) & (prp_mssng[-1*length(prp_mssng)] < ctff)) + 1
  ends <- which((prp_mssng[-1] < ctff) & (prp_mssng[-1*length(prp_mssng)] > ctff))
  lngths <- ends - strts
  pull <- lngths > mn_lngth
  strts <- strts[pull]
  ends <- ends[pull]
  new_sqs <- vector('list', length=length(strts))
  for(i in seq_along(strts)) {
    strt <- strts[i]
    end <- ends[i]
    new_sq <- vector('list', length=length(sq))
    for(j in seq_along(sq)) {
      new_sq[[j]] <- substr(sq[[j]], start=strt, stop=end)
    }
    names(new_sq) <- names(sq)
    new_sqs[[i]] <- new_sq
  }
  new_sqs
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