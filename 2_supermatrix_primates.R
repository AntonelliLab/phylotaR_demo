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

# CHECK AND WRITE OUT
all(sapply(sprmtrx, nchar) == nchar(sprmtrx[[1]]))
writeSqs(sprmtrx, fl=file.path(wd, 'supermatrix.fasta'))

# RAxML
# Warning: partition.txt may need minor modification depending on gene type
inpt <- file.path(wd, 'supermatrix.fasta')
system(paste0('raxmlHPC -m GTRGAMMA -T 2 -f a -N 10 -p 1234 -x 1234 -n primates -s ', inpt))

# CLEAN-UP
file.rename('RAxML_bestTree.primates', file.path(wd, 'tree.tre'))
file.remove(list.files(pattern='RAxML*'))
