library(phylotaR)
wd <- 'primates'
tdpth <- 'taxdump.tar.gz'
ncbi_dr <- file.path('ncbi-blast-2.7.1+', 'bin')
txid <- 9443  # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9443
setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr, tdpth=tdpth, ncps=1, mxsqs=10000)
run(wd=wd, nstages=4)
