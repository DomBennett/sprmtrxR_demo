# Libs ----
library(phylotaR)

# Vars ----
wd <- 'primates'
ncbi_dr <- '/usr/bin/'
txid <- 376911  # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9443
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
run(wd = wd)
