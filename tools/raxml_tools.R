seqs_write <- function(seqs, fl) {
  fasta <- ''
  for (i in seq_along(seqs)) {
    sq <- seqs[[i]]
    dfln <- names(seqs)[[i]]
    dfln <- gsub(pattern = '\\s', replacement = '_', x = dfln)
    fasta <- paste0(fasta, '>', dfln, '\n', sq, '\n\n')
  }
  cat(fasta, file = fl)
}

partition <- function(lngths, fl) {
  gene <- strt <- 1
  prttn_txt <- ''
  for (lngth in lngths) {
    end <- lngth + strt - 1
    prttn_txt <- paste0(prttn_txt, 'DNA, gene', gene, ' = ', strt, '-', end,
                        '\n')
    strt <- end + 1
    gene <- gene + 1
  }
  cat(prttn_txt, file = fl)
}
