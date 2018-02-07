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