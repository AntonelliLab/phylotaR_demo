# Compare outputs from PhyLoTa and phylotaR
# as per reviewer suggestions
# LIBS
library(phylotaR)

# PALMS ----
all_cls <- read_phylota('palms')
# get phylotar_gis
cids <- all_cls@cids
phylotar_gis <- vector(mode = 'list', length = length(cids))
names(phylotar_gis) <- cids
for (cid in cids) {
  phylotar_gis[[cid]] <- get_sq_slot(phylota = all_cls, cid = cid, slt_nm = 'gi')
}
# get phylota_gis
pgfile <- file.path('results', 'palms_phylota_gis.RData')
if (file.exists(pgfile)) {
  phylota_gis <- readRDS(pgfile)
} else {
  # http://phylota.net/cgi-bin/sql_getdesc.cgi?c=0&ti=169703&mode=0&db=194
  phylota_cids <- as.character(0:67)
  phylota_gis <- vector(mode = 'list', length = length(phylota_cids))
  names(phylota_gis) <- phylota_cids
  for (cid in phylota_cids) {
    cat(cid, '\n')
    txid <- 169703
    url <- paste0('http://phylota.net/cgi-bin/sql_getcluster_fasta.cgi?format=gi&db=194&ti=',
                  txid,'&cl=', cid, '&ntype=1')
    fasta <- readLines(url, warn = FALSE)
    gi_lines <- fasta[grepl(pattern = '>gi', x = fasta)]
    phylota_gis[[cid]] <- gsub(pattern = '[^0-9]+', replacement = '', x = gi_lines)
  }
  saveRDS(phylota_gis, file = file.path('results', 'palms_phylota_gis.RData'))
}
# compare
res <- matrix(data = NA, nrow = length(phylota_gis), ncol = length(phylotar_gis))
for (i in seq_along(phylota_gis)) {
  gis1 <- phylota_gis[[i]]
  for (j in seq_along(phylotar_gis)) {
    gis2 <- phylotar_gis[[j]]
    res[i, j] <- sum(gis1 %in% gis2)/length(gis1)
  }
}
max_prp <- apply(X = res, MARGIN = 1, max)
mean(max_prp[sapply(phylota_gis, length) >= 4])
# 0.9202058 of PhyLoTa cluster sequences are in phylotaR
slctd <- c('0', '2', '4', '16', '1', '10', '3', '35', '93', '139')
slctd_res <- res[ , names(phylotar_gis) %in% slctd]
mean(apply(X = slctd_res, MARGIN = 2, max))
# 0.9919231 for selected clusters

# PRIMATES ----
all_cls <- read_phylota('primates')
n_taxa <- get_cl_slot(all_cls, cid = all_cls@cids, slt_nm = 'ntx')
cltyps <- get_cl_slot(all_cls, all_cls@cids, slt_nm = 'typ')
picls <- n_taxa >= 4
picls <- all_cls@cids[picls]
# get phylotar_gis ====
phylotar_gis <- vector(mode = 'list', length = length(picls))
names(phylotar_gis) <- picls
for (i in seq_along(picls)) {
  cat(i, '/', length(picls), '\n', sep = '')
  phylotar_gis[[picls[[i]]]] <- get_sq_slot(phylota = all_cls, cid = picls[[i]],
                                            slt_nm = 'gi')
}
# get phylotar_tis ====
phylotar_tis <- vector(mode = 'list', length = length(picls))
names(phylotar_tis) <- picls
for (i in seq_along(picls)) {
  cat(i, '/', length(picls), '\n', sep = '')
  phylotar_tis[[picls[[i]]]] <- get_sq_slot(phylota = all_cls, cid = picls[[i]],
                                            slt_nm = 'txid')
}
# get phylota_gis ====
pgfile <- file.path('results', 'primates_phylota_gis.RData')
if (file.exists(pgfile)) {
  phylota_gis <- readRDS(pgfile)
} else {
  # http://phylota.net/cgi-bin/sql_getdesc.cgi?c=0&ti=9443&mode=0&db=194
  phylota_cids <- as.character(0:2022)
  phylota_gis <- vector(mode = 'list', length = length(phylota_cids))
  names(phylota_gis) <- phylota_cids
  for (cid in phylota_cids) {
    cat(cid, '\n')
    txid <- 9443
    url <- paste0('http://phylota.net/cgi-bin/sql_getcluster_fasta.cgi?format=gi&db=194&ti=',
                  txid,'&cl=', cid, '&ntype=1')
    fasta <- readLines(url, warn = FALSE)
    gi_lines <- fasta[grepl(pattern = '>gi', x = fasta)]
    phylota_gis[[cid]] <- gsub(pattern = '[^0-9]+', replacement = '', x = gi_lines)
  }
  saveRDS(phylota_gis, file = pgfile)
}
# get phylota_tis ====
ptfile <- file.path('results', 'primates_phylota_tis.RData')
if (file.exists(ptfile)) {
  phylota_tis <- readRDS(ptfile)
} else {
  # http://phylota.net/cgi-bin/sql_getdesc.cgi?c=0&ti=9443&mode=0&db=194
  phylota_cids <- as.character(0:2022)
  phylota_tis <- vector(mode = 'list', length = length(phylota_cids))
  names(phylota_tis) <- phylota_cids
  for (cid in phylota_cids) {
    cat(cid, '\n')
    txid <- 9443
    url <- paste0('http://phylota.net/cgi-bin/sql_getcluster_fasta.cgi?format=ti&db=194&ti=',
                  txid,'&cl=', cid, '&ntype=1')
    fasta <- readLines(url, warn = FALSE)
    ti_lines <- fasta[grepl(pattern = '>ti', x = fasta)]
    phylota_tis[[cid]] <- gsub(pattern = '[^0-9]+', replacement = '', x = ti_lines)
  }
  saveRDS(phylota_tis, file = ptfile)
}
# get phylota_lngths ====
plfile <- file.path('results', 'primates_phylota_lngths.RData')
if (file.exists(plfile)) {
  phylota_lngths <- readRDS(plfile)
} else {
  # http://phylota.net/cgi-bin/sql_getdesc.cgi?c=0&ti=9443&mode=0&db=194
  phylota_cids <- as.character(0:2022)
  phylota_lngths <- vector(mode = 'list', length = length(phylota_cids))
  names(phylota_lngths) <- phylota_cids
  for (cid in phylota_cids) {
    cat(cid, '\n')
    txid <- 9443
    url <- paste0('http://phylota.net/cgi-bin/sql_getcluster_fasta.cgi?format=gi&db=194&ti=',
                  txid,'&cl=', cid, '&ntype=1')
    fasta <- readLines(url, warn = FALSE)
    gilines <- c(which(grepl(pattern = '>gi', x = fasta)), length(fasta))
    res <- list()
    for (i in 1:(length(gilines) - 1)) {
      gi <- gsub(pattern = '[^0-9]+', replacement = '', x = fasta[[gilines[[i]]]])
      sqis <- (gilines[[i]] + 1):(gilines[[i + 1]] - 1)
      res[[gi]] <- sum(sapply(fasta[sqis], nchar))
    }
    phylota_lngths[[cid]] <- res
  }
  saveRDS(phylota_lngths, file = plfile)
}
# work out model organisms ====
all_tis <- unique(unlist(phylotar_tis))
ps <- phylotaR:::ldPrmtrs('primates')
ps[['wd']] <- file.path(getwd(), 'primates')
ps[['v']] <- TRUE
nsqs_txids <- rep(NA, length(all_tis))
names(nsqs_txids) <- all_tis
for (i in seq_along(all_tis)) {
  cat(i, '/', length(all_tis), '\n', sep = '')
  nsqs_txids[[all_tis[[i]]]] <- phylotaR:::nSqs(txid = all_tis[[i]],
                                                ps = ps)
}
mdl_orgs <- names(nsqs_txids)[nsqs_txids > parameters()[['mdlthrs']]]
# compare ====
# without model organisms, controlling for differences in taxonomy and seq length
res <- matrix(data = NA, nrow = length(phylota_gis), ncol = length(phylotar_gis))
for (i in seq_along(phylota_gis)) {
  gis1 <- phylota_gis[[i]]
  tis1 <- phylota_tis[[i]]
  lngths <- as.numeric(phylota_lngths[[i]])
  pulllngths <- lngths > 250 & lngths < 2000
  gis1 <- gis1[pulllngths]
  tis1 <- tis1[pulllngths]
  pullmo <- !tis1 %in% mdl_orgs
  gis1 <- gis1[pullmo]
  tis1 <- tis1[pullmo]
  for (j in seq_along(phylotar_gis)) {
    gis2 <- phylotar_gis[[j]]
    tis2 <- phylotar_tis[[j]]
    pull1 <- tis1 %in% tis2
    pull2 <- tis2 %in% tis1
    res[i, j] <- sum(gis1[pull1] %in% gis2[pull2])/sum(pull1)
  }
}
res[is.na(res)] <- 0
max_prp <- apply(X = res, MARGIN = 1, max)
mean(max_prp[sapply(phylota_gis, length) >= 4])
# 0.5512457 of PhyLoTa cluster sequences are in phylotaR
slctd <- c('1', '8', '4', '13', '12', '30', '28', '37', '40', '43')
slctd_res <- res[ , names(phylotar_gis) %in% slctd]
mean(apply(X = slctd_res, MARGIN = 2, max))
# 0.9996947 for selected clusters

# get missing ====
# get all gis not in phylotaR results
all_tis <- unique(unlist(phylotar_tis))
all_gis <- NULL
for (i in seq_along(phylota_gis)) {
  gis <- phylota_gis[[i]]
  tis <- phylota_tis[[i]]
  lngths <- as.numeric(phylota_lngths[[i]])
  pulllngths <- lngths > 250 & lngths < 2000
  gis <- gis[pulllngths]
  tis <- tis[pulllngths]
  pullmo <- !tis %in% mdl_orgs
  gis <- gis[pullmo]
  tis <- tis[pullmo]
  all_gis <- c(gis[tis %in% all_tis], all_gis)
}
all_gis <- unique(all_gis)
mssng <- all_gis[!all_gis %in% unique(unlist(phylotar_gis))]
mssng_rcrds <- rentrez::entrez_summary(db = 'nucleotide', id = mssng[1:300])
btchs <- seq(301, length(mssng), 300)
for (i in seq_along(btchs)[-1]) {
  cat(i, '/', length(btchs), '\n', sep = '')
  btch <- mssng[btchs[(i - 1)]:btchs[i]]
  mssng_rcrds <- c(mssng_rcrds, rentrez::entrez_summary(db = 'nucleotide', id = btch))
}
length(mssng_rcrds)
length(mssng)
defunct <- sapply(mssng_rcrds, '[[', 'status') == 'replaced'
sum(defunct)/length(mssng)
# 0.09630419 are defunct IDs
accessions <- sapply(mssng_rcrds, '[[', 'accessionversion')
accessions <- sub('\\.[0-9]$', '', accessions)
in_small_clusters <- vapply(X = accessions, FUN = function(x) {
  pttrn <- paste0('^', x, '.*')
  any(grepl(pattern = pttrn, x = all_cls@sids))
}, logical(1))
sum(in_small_clusters)/length(mssng)
example_small <- mssng_rcrds[[mssng[in_small_clusters][[1]]]][['accessionversion']]
for (cid in all_cls@cids) {
  cl <- all_cls@cls[[cid]]
  if (example_small %in% cl@sids) {
    stop()
  }
}
# 0.2999289 are in clusters that are not phylogenetically informative
sqpth <- file.path('primates', 'cache', 'sqs')
sqfls <- list.files(path = sqpth)
all_sqs <- NULL
for (sqfl in sqfls) {
  sqs <- readRDS(file = file.path(sqpth, sqfl))
  all_sqs <- c(all_sqs, sqs@ids)
}
downloaded <- vapply(X = accessions, FUN = function(x) {
  pttrn <- paste0('^', x, '.*')
  any(grepl(pattern = pttrn, x = all_sqs))
}, logical(1))
sum(!downloaded & !defunct)/length(mssng)
# 0.002487562 were not dowloaded
unclustered <- mssng[downloaded & !defunct & !in_small_clusters]
length(unclustered)/length(mssng)
# 0.5707178 downloaded but not in any clusters
str(mssng_rcrds[[unclustered[[1]]]])
unclustered_lengths <- unlist(sapply(mssng_rcrds[unclustered], '[[', 'slen'))
median(unclustered_lengths)
clustered_lengths <- unlist(sapply(mssng_rcrds[in_small_clusters], '[[', 'slen'))
median(clustered_lengths)


# i <- 1
# j <- 7 #10
# gis1 <- phylota_gis[[i]]
# tis1 <- phylota_tis[[i]]
# pull <- !tis1 %in% c('9606', '9600', '9544', '9598', '9601')
# gis1 <- gis1[pull]
# tis1 <- tis1[pull]
# gis2 <- phylotar_gis[[j]]
# sum(gis1 %in% gis2)/length(gis1)
# mssng7 <- gis1[!gis1 %in% gis2]
# 
# 
# mssng10 %in% mssng7
# '32435531'
# which(grepl(pattern = '^DQ495311.*', x = all_cls@sids))
# all_cls@txdct@rcrds[['9551']]
# cl <- all_cls@cls@cls[[picls[10]]]
# '9551' %in% cl@txids
# '1411634' %in% cl@txids
# 
# sqpth <- file.path('primates', 'cache', 'sqs')
# sqfls <- list.files(path = sqpth)
# for (sqfl in sqfls) {
#   sqs <- readRDS(file = file.path(sqpth, sqfl))
#   if (any(grepl(pattern = '^AY665629.*', x = sqs@ids))) {
#     stop()
#   }
#   # if ('1411634' %in% sqs@txids) {
#   #   stop()
#   # }
# }
