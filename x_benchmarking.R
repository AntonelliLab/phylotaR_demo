# Demonstrating how different parameters impact the results
library(phylotaR)

# FUNCTIONS ----
dir_setup <- function(wd) {
  if (dir.exists(wd)) {
    unlink(wd, recursive = TRUE)
  }
  dir.create(wd)
}
log_times_get <- function(fldr, base_wd) {
  lgfl <- file.path(base_wd, fldr, 'log.txt')
  lines <- readLines(con = lgfl)
  stage_tms <- stage_ends <- stage_starts <- c('taxise' = NA,
                                               'download' = NA,
                                               'cluster' = NA,
                                               'cluster\\^2' = NA)
  stage_nms <- names(stage_tms)
  for (ln in lines) {
    for (stgnm in stage_nms) {
      pttrn <- paste0('Starting stage ', stgnm, ': ')
      if (grepl(pattern = pttrn, x = ln, ignore.case = TRUE)) {
        ln <- sub(pattern = pttrn,
                  replacement = '', x = ln, ignore.case = TRUE)
        ln <- gsub(pattern = '(\\[|\\])', replacement = '', x = ln)
        stage_starts[[stgnm]] <- ln
      }
      pttrn <- paste0('Completed stage ', stgnm, ': ')
      if (grepl(pttrn, ln, ignore.case = TRUE)) {
        ln <- sub(pattern = pttrn,
                  replacement = '', x = ln, ignore.case = TRUE)
        ln <- gsub(pattern = '(\\[|\\])', replacement = '', x = ln)
        stage_ends[[stgnm]] <- ln
      }
    }
  }
  for (stgnm in stage_nms) {
    stage_tms[[stgnm]] <- difftime(as.POSIXct(stage_ends[[stgnm]]),
                                   as.POSIXct(stage_starts[[stgnm]]),
                                   units = 'mins')
  }
  stage_tms
}

# RUN ----
# shared vars
base_wd <- 'benchmarking'
txid <- 4710  # Areaceae
ncbi_dr <- file.path('ncbi-blast-2.7.1+', 'bin')

cat('.... defaults\n')
wd <- file.path(base_wd, 'defaults')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
run(wd = wd)

cat('.... mnsql_100\n')
wd <- file.path(base_wd, 'mnsql_100')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mnsql = 100)
run(wd = wd)

cat('.... mnsql_500\n')
wd <- file.path(base_wd, 'mnsql_500')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mnsql = 500)
run(wd = wd)

cat('.... mxsql_1000\n')
wd <- file.path(base_wd, 'mxsql_1000')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mxsql = 1000)
run(wd = wd)

cat('.... mxsql_5000\n')
wd <- file.path(base_wd, 'mxsql_5000')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mxsql = 5000)
run(wd = wd)

cat('.... mnsql_500_mxsql_1000\n')
wd <- file.path(base_wd, 'mnsql_500_mxsql_1000')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mnsql = 500, mxsql = 1000)
run(wd = wd)

cat('.... mnsql_100_mxsql_5000\n')
wd <- file.path(base_wd, 'mnsql_100_mxsql_5000')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mnsql = 100, mxsql = 5000)
run(wd = wd)

cat('.... mxsqs_500\n')
wd <- file.path(base_wd, 'mxsqs_500')
dir_setup(wd)
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr, mxsqs = 500)
run(wd = wd)

# STATS ----
fldrs <- list.files(base_wd)
fldrs <- fldrs[-length(fldrs)]
stat_nms <- c('tm', 'taxise_tm', 'download_tm', 'cluster_tm',
              'cluster2_tm', 'nsqs', 'ncls', 'nmerged', 'nsubtree',
              'ndirect', 'nparaphyly', 'ntx', 'median_tx_per_cl',
              'median_sq_per_cl', 'median_sqlngth')
stats <- matrix(NA, nrow = length(fldrs), ncol = length(stat_nms))
colnames(stats) <- stat_nms
rownames(stats) <- fldrs
for (fldr in fldrs) {
  cls <- read_phylota(file.path(base_wd, fldr))
  log_times <- log_times_get(fldr = fldr, base_wd = base_wd)
  stats[fldr, 'tm'] <- sum(log_times)
  stats[fldr, 'taxise_tm'] <- log_times[['taxise']]
  stats[fldr, 'download_tm'] <- log_times[['download']]
  stats[fldr, 'cluster_tm'] <- log_times[['cluster']]
  stats[fldr, 'cluster2_tm'] <- log_times[['cluster\\^2']]
  stats[fldr, 'nsqs'] <- length(cls@sids)
  stats[fldr, 'ncls'] <- length(cls@cids)
  cltyps <- get_cl_slot(cls, cls@cids, slt_nm = 'typ')
  stats[fldr, 'nmerged'] <- sum(cltyps == 'merged')
  stats[fldr, 'nsubtree'] <- sum(cltyps == 'subtree')
  stats[fldr, 'ndirect'] <- sum(cltyps == 'direct')
  stats[fldr, 'nparaphyly'] <- sum(cltyps == 'paraphyly')
  stats[fldr, 'ntx'] <- length(unique(cls@txids))
  stats[fldr, 'median_tx_per_cl'] <-
    median(vapply(cls@cids, function(cid) {
      length(unique(cls@cls[[cid]]@txids))},
      integer(1)))
  stats[fldr, 'median_sq_per_cl'] <-
    median(vapply(cls@cids, function(cid) {
      length(cls@cls[[cid]]@sids)},
      integer(1)))
  stats[fldr, 'median_sqlngth'] <-
    median(vapply(cls@sids, function(sid) {
      cls@sqs[[sid]]@nncltds},
      integer(1)))
}

# OUTPUT
write.csv(x = stats, file = file.path('results', 'benchmarking.csv'))
