# Combine alignments into a supermatrix, run RAxML

# LIBS
source('tools', 'supermatrix_tools.R')

# VARS
wd <- 'palms'

# INPUT
sqs1 <- readSqs(file.path(wd, 'alignment1.fasta'))
sqs2 <- readSqs(file.path(wd, 'alignment2.fasta'))

# GEN PARTITION TEXT
lngths <- c(nchar(sqs1[[1]]), nchar(sqs2[[1]]))
partition(lngths, fl=file.path(wd, 'partition.txt'))

# GEN SUPERMARTIX
fllr1 <- paste0(rep('-', lngths[[1]]), collapse='')
fllr2 <- paste0(rep('-', lngths[[2]]), collapse='')
all_nms <- unique(c(names(sqs1), names(sqs2)))
all_nms <- sort(all_nms)
sprmtrx <- vector('list', length=length(all_nms))
names(sprmtrx) <- all_nms
for(nm in all_nms) {
  sq1 <- sqs1[[nm]]
  sq1 <- ifelse(is.null(sq1), fllr1, sq1)
  sq2 <- sqs2[[nm]]
  sq2 <- ifelse(is.null(sq2), fllr2, sq2)
  sprmtrx[[nm]] <- paste0(sq1, sq2)
}

# CHECK AND WRITE OUT
all(sapply(sprmtrx, nchar) == nchar(sprmtrx[[1]]))
sprmtrx <- sprmtrx[!grepl('^Arecaceae', names(sprmtrx))]
writeSqs(sprmtrx, fl=file.path(wd, 'supermatrix.fasta'))

# RAxML
# Warning: partition.txt may need minor modification depending on gene type
inpt <- file.path(wd, 'supermatrix.fasta')
system(paste0('raxmlHPC -m GTRGAMMA -f a -N 10 -p 1234 -x 1234 -n palms -s ', inpt))

# CLEAN-UP
file.rename('RAxML_bestTree.palms', file.path(wd, 'tree.tre'))
file.remove(list.files(pattern='RAxML*'))