# Generate supermatrices for mono groups
# Generate combined super-supermatrix

# Load ----
source(file = file.path('tools', 'load.R'))

# Vars ----
column_cutoff <- .51
tip_cutoff <- .75
alignment_dir <- file.path('primates', 'alignments')
supermatrix_dir <- file.path('primates', 'supermatrix')
if (!dir.exists(supermatrix_dir)) {
  dir.create(supermatrix_dir)
}
groups_dir <- file.path('primates', 'groups')

# Data ----
groups <- readRDS(file = file.path(groups_dir, 'groups.RData'))
alfls <- file.path(alignment_dir, list.files(path = alignment_dir,
                                             pattern = '.fasta'))

# Alignments ----
alignments <- list()
for (alfl in alfls) {
  alignments[[alfl]] <- seqs_read(fl = alfl)
}
names(alignments) <- gsub(pattern = '[^0-9]', replacement = '',
                          x = names(alignments))

# Supermatrices ----
smatrices <- smatrices_get(groups = groups, alignments = alignments,
                           column_cutoff = column_cutoff,
                           tip_cutoff = tip_cutoff)

# Write out smatrices ----
for (i in seq_along(smatrices)) {
  smatrix <- smatrices[[i]]
  grp_id <- names(smatrices)[[i]]
  seqs_write(seqs = smatrix, fl = file.path(supermatrix_dir,
                                            paste0(grp_id, '.fasta')))
  partition_file <- file.path(supermatrix_dir, paste0(grp_id, '_partition.txt'))
  partition(lngths = attr(smatrix, 'nbps'), fl = partition_file)
}
