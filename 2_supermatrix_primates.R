# Combine alignments into a supermatrix, run RAxML

# LIBS
source(file.path('tools', 'supermatrix_tools.R'))

# VARS
wd <- 'primates'
prosimians <- read.delim(file.path('expected', 'prosimians_list.txt'),
                         stringsAsFactors=FALSE)[ ,1]

# INPUT
alfls <- list.files(wd, pattern='alignment[0-9]+.fasta')
als <- vector('list', length=length(alfls))
for(i in seq_along(alfls)) {
  als[[i]] <- readSqs(file.path(wd, alfls[i]))
}

# DROP OVERHANGING EDGES
(sapply(als, function(x) nchar(x[[1]])))
als <- drpOvrhngs(als)
(sapply(als, function(x) nchar(x[[1]])))
n_taxa <- sapply(als, length)

# GEN PARTITION TEXT
lngths <- sapply(als, function(x) nchar(x[[1]]))
partition(lngths, fl=file.path(wd, 'partition.txt'))

# GEN SUPERMARTIX
fllrs <- sapply(lngths, function(x) paste0(rep('-', x),
                                           collapse=''))
all_nms <- unique(unlist(sapply(als, names)))
all_nms <- sort(all_nms)
pull <- !grepl('\\ssp\\.', all_nms)
all_nms <- all_nms[pull]
sprmtrx <- vector('list', length=length(all_nms))
names(sprmtrx) <- all_nms
for(nm in all_nms) {
  al <- ''
  for(i in seq_along(als)) {
    tmp <- als[[i]][[nm]]
    tmp <- ifelse(is.null(tmp), fllrs[[i]], tmp)
    al <- paste0(al, tmp)
  }
  sprmtrx[[nm]] <- al
}

# DROP TAXA WITH TOO MANY GAPS
ngaps <- sapply(gregexpr('-', sprmtrx), length)
pull <- ngaps < nchar(sprmtrx[[1]])
sprmtrx <- sprmtrx[pull]

# CHECK AND WRITE OUT
all(sapply(sprmtrx, nchar) == nchar(sprmtrx[[1]]))
names(sprmtrx) <- gsub('\\s', '_', names(sprmtrx))
writeSqs(sprmtrx, fl=file.path(wd, 'supermatrix.fasta'))

# OUTGROUP
nms <- sub('_[0-9]', '', names(sprmtrx))
pull <- sapply(nms, function(x) any(grepl(x, prosimians)))
outgroup <- names(sprmtrx)[pull]
outgroup <- paste0(outgroup, collapse=',')

# RAxML
# Warning: partition.txt may need minor modification depending on gene type
inpt <- file.path(wd, 'supermatrix.fasta')
prttnfl <- file.path(wd, 'partition.txt')
system(paste0('raxmlHPC -f a -m GTRGAMMA -T 2 -# 10 -p ',
              sample(0:10000000, 1), ' -x ', sample(0:10000000, 1),
              ' -n primates -s ', inpt, ' -q ', prttnfl))#, ' -o ', outgroup))
# consensus
system('raxmlHPC -m GTRCAT -J MR -z RAxML_bootstrap.primates -n primates_con')

# CLEAN-UP
file.rename('RAxML_bestTree.primates', file.path(wd, 'best_tree.tre'))
file.rename('RAxML_bootstrap.primates', file.path(wd, 'bootstraps.tre'))
file.rename('RAxML_MajorityRuleConsensusTree.primates_con',
            file.path(wd, 'consensus.tre'))
file.remove(list.files(pattern='RAxML*'))
