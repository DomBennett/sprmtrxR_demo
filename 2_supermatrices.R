# Generate supermatrices for mono groups
# Generate combined super-supermatrix

# Load ----
source(file = file.path('tools', 'load.R'))

# Vars ----
alignment_dir <- file.path('primates', 'alignments')
groups_dir <- file.path('primates', 'groups')
groups <- readRDS(file = file.path(groups_dir, 'groups.RData'))
alfls <- file.path(alignment_dir, list.files(path = alignment_dir,
                                             pattern = '.fasta'))

alignments <- list()
for (alfl in alfls) {
  alignments[[alfl]] <- seqs_read(fl = alfl)
}

smatrices_get <- function(groups, alignments, column_cutoff = .75,
                          tip_cutoff = column_cutoff, min_ntips = 5,
                          min_ngenes = 5) {
  nbps <- vapply(X = alignments, FUN = function(x) length(x[[1]]),
                 FUN.VALUE = integer(1))
  smatrices <- vector(mode = 'list', length = length(groups))
  names(smatrices) <- names(groups)
  for (grp_id in names(groups)) {
    nms <- groups[[grp_id]]
    smatrix <- smatrix_get(nms = nms, alignments = alignments, nbps = nbps)
    smatrix <- drop_columns(smatrix = smatrix, cutoff = column_cutoff)
    smatrix <- drop_tips(smatrix = smatrix, cutoff = tip_cutoff)
    if (smatrix) {
      smatrices[[grp_id]] <- smatrix
    }
  }
  res <- list('smatrices' = smatrices, 'groups' = names(groups))
  class(res) <- 'smatrices'
  res
}
