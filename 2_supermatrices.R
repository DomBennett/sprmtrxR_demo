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
names(alignments) <- gsub(pattern = '[^0-9]', replacement = '',
                          x = names(alignments))

smatrices_get <- function(groups, alignments, column_cutoff = .75,
                          tip_cutoff = column_cutoff, min_ntips = 5,
                          min_ngenes = 5, min_nbps = 200) {
  res <- list()
  all_tips <- character(0)
  for (grp_id in names(groups)) {
    nms <- groups[[grp_id]]
    # select sequences from alignments
    seqs_list <- seqs_select(nms = nms, seqs_list = alignments)
    # filter selected sequences
    seqs_list <- seqs_filter(seqs_list = seqs_list, cutoff = column_cutoff,
                             min_nbps = min_nbps)
    # merge into supermatrix
    smatrix <- smatrix_get(seqs_list = seqs_list)
    # drop tips
    smatrix <- drop_tips(smatrix = smatrix, cutoff = tip_cutoff)
    if (length(smatrix) >= min_ntips &
        length(attr(smatrix, 'genes')) >= min_ngenes) {
      res[[grp_id]] <- smatrix
      # record tips
      all_tips <- c(all_tips, names(smatrix))
    }
  }
  attr(res, 'tips') <- all_tips
  class(res) <- 'res'
  res
}

smatrices <- res
i <- 1

# Write out smatrices ----
for (i in seq_along(smatrices)) {
  smatrix <- smatrices[[i]]
  # different lengths
  smatrix <- smatrix[-1*c(8,9)]
  grp_id <- names(smatrices)[[i]]
  seqs_write(seqs = smatrix, fl = paste0(grp_id, '.fasta'))
  partition(lngths = attr(smatrix, 'nbps'), fl = 'partition.txt')
}
