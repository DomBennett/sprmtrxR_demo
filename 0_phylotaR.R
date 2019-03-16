# Libs ----
library(phylotaR)
library(outsider)

# Mafft import
# module_install('dombennett/om..mafft')
mafft <- module_import('mafft', 'dombennett/om..mafft')

# Vars ----
wd <- 'primates'
ncbi_dr <- '/usr/local/ncbi/blast/bin/'
txid <- 376911  # https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9443

# Run ----
setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
run(wd = wd)

# Parse ----
primates <- read_phylota(wd)
primates <- drop_by_rank(primates, rnk = 'species', n = 1,
                         choose_by = c("pambgs", "age", "nncltds"),
                         greatest = c(FALSE, FALSE, TRUE))
ntaxa <- get_ntaxa(phylota = primates, cid = primates@cids)
cids <- names(ntaxa)[ntaxa > 10]

# Output
sqfl_wd <- file.path(wd, 'clusters')
dir.create(sqfl_wd)
for (i in seq_along(cids)) {
  cid <- cids[i]
  sids <- primates@clstrs[[cid]]@sids
  txids <- get_txids(primates, cid = cid, rnk = 'species')
  scnms <- get_tx_slot(primates, txids, 'scnm')
  n <- sapply(seq_along(scnms), function(x) 
    sum(scnms[x] == scnms[x:length(scnms)]))
  infile <- file.path(sqfl_wd, paste0('sequences', i, '.fasta'))
  write_sqs(phylota = primates, outfile = infile, sid = sids, sq_nm = scnms)
}

# Alignment
sqfls <- list.files(sqfl_wd)
for (sqfl in sqfls) {
  alfl <- sub(pattern = 'sequences', replacement = 'alignment', x = sqfl)
  mafft('--auto', file.path(sqfl_wd, sqfl), '>', file.path(sqfl_wd, alfl))
}
