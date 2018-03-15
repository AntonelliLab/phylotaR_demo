library(phylotaR)
wd <- 'palms'
ncbi_dr <- file.path('ncbi-blast-2.7.1+', 'bin')
txid <- 4710  # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4710
setUp(wd=wd, txid=txid, ncbi_dr=ncbi_dr)
run(wd=wd)
