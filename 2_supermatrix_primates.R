# Combine alignments into a supermatrix, run RAxML

# LIBS
source(file.path('tools', 'supermatrix_tools.R'))

# VARS
wd <- 'primates'

# INPUT
sqfls <- list.files(wd, pattern='alignment[0-9]+.fasta')
sqs <- vector('list', length=length(sqfls))
for(i in seq_along(sqfls)) {
  sqs[[i]] <- readSqs(file.path(wd, sqfls[i]))
}

# DROP OVERHANGING EDGES
sqs <- drpOvrhngs(sqs)
(sapply(sqs, function(x) nchar(x[[1]])))
# drop first cluster -- too big
new_sqs <- spltUp(sqs[[1]], mn_lngth=700)
(sapply(new_sqs, function(x) nchar(x[[1]])))
sqs <- sqs[-1]
sqs <- c(sqs, new_sqs)

# GEN PARTITION TEXT
lngths <- sapply(sqs, function(x) nchar(x[[1]]))
partition(lngths, fl=file.path(wd, 'partition.txt'))

# GEN SUPERMARTIX
fllrs <- sapply(lngths, function(x) paste0(rep('-', x),
                                           collapse=''))
all_nms <- unique(unlist(sapply(sqs, names)))
all_nms <- sort(all_nms)
sprmtrx <- vector('list', length=length(all_nms))
names(sprmtrx) <- all_nms
for(nm in all_nms) {
  sq <- ''
  for(i in seq_along(sqs)) {
    tmp <- sqs[[i]][[nm]]
    tmp <- ifelse(is.null(tmp), fllrs[[i]], tmp)
    sq <- paste0(sq, tmp)
  }
  sprmtrx[[nm]] <- sq
}

# DROP TAXA WITH TOO MANY GAPS
ngaps <- sapply(gregexpr('-', sprmtrx), length)
pull <- ngaps < nchar(sprmtrx[[1]])
sprmtrx <- sprmtrx[pull]

# CHECK AND WRITE OUT
all(sapply(sprmtrx, nchar) == nchar(sprmtrx[[1]]))
writeSqs(sprmtrx, fl=file.path(wd, 'supermatrix.fasta'))

# DEFINE OUTGROUP
prosimians <- c('Hapalemur', 'Lemur', 'Varecia', 'Eulemur',
                'Avahi', 'Propithecus', 'Mirza', 'Daubentonia',
                'Galago', 'Arctocebus', 'Perodicticus', 'Loris',
                'Nycticebus', 'Galagoides', 'Euoticus')
outgroup <- prosimians[prosimians %in% names(sprmtrx)]
outgroup <- paste0(outgroup, collapse=',')

# RAxML
# Warning: partition.txt may need minor modification depending on gene type
inpt <- file.path(wd, 'supermatrix.fasta')
system(paste0('raxmlHPC -f a -m GTRGAMMA -T 2 -# 100 -p ',
              sample(0:10000000, 1), ' -x ', sample(0:10000000, 1),
              ' -n primates -s ', inpt))
# consensus
system('raxmlHPC -m GTRCAT -J MR -z RAxML_bootstrap.primates -n primates_con')

# CLEAN-UP
file.rename('RAxML_bestTree.primates', file.path(wd, 'best_tree.tre'))
file.rename('RAxML_bootstrap.primates', file.path(wd, 'bootstraps.tre'))
file.rename('RAxML_MajorityRuleConsensusTree.primates_con',
            file.path(wd, 'consensus.tre'))
file.remove(list.files(pattern='RAxML*'))
