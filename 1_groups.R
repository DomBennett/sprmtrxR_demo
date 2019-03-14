# Load ----
source(file = file.path('tools', 'load.R'))

# Vars ----
max_size <- 60
min_size <- 5
input_dir <- file.path('primates', 'alignments')
output_dir <- file.path('primates', 'groups')
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
alfls <- file.path(input_dir, list.files(path = input_dir, pattern = '.fasta'))

# Data ----
tree <- treeman::readTree(file = file.path('0_data', 'primate.tre'))
tree_nms <- tree@tips
tree_pttrns <- sub(pattern = '_.*$', replacement = '', x = tree_nms)
algnmnt_nms <- nms_from_alignments(flpths = alfls)
algnmnt_pttrns <- sub(pattern = '\\s.*$', replacement = '', x = algnmnt_nms)
matched_nms <- nm_match(algnmnt_nms = algnmnt_nms, tree_nms = tree_nms,
                        algnmnt_pttrns = algnmnt_pttrns,
                        tree_pttrns = tree_pttrns, max_dist = .15)
groups <- groups_get(matched_nms = matched_nms, tree = tree,
                     max_size = max_size, min_size = min_size)

# Output ----
output_file <- file.path(output_dir, 'groups.RData')
saveRDS(object = groups, file = output_file)

