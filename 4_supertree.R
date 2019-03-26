# take subtrees and add to supertree

# Vars
input_dir <- file.path('primates', 'raxml')
tree_files <- file.path(input_dir, list.files(input_dir, 'bestTree'))


ids <- unlist(lapply(X = strsplit(x = tree_files, split = '\\.'),
                     FUN = '[[', 2))
names(tree_files) <- ids

trstrs <- vapply(X = tree_files, FUN = readLines, FUN.VALUE = character(1))

strsplit(x = trstrs[['unmatched']], split = ',')

res <- trstrs[['unmatched']]
ids <- ids[ids != 'unmatched']
for (id in ids) {
  matching <- grepl(pattern = id, x = res)
  if (matching) {
    split_res <- strsplit(x = res, split = id)[[1]]
    unrooted <- sub(pattern = ':[0-9\\.]+;', replacement = '', x = trstrs[[id]])
    res <- paste0(split_res[[1]], unrooted, split_res[[2]])
  }
}

tree <- ape::read.tree(text = res)
plot(tree)
