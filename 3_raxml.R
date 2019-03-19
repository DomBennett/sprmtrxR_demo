
# TODO:
# - ensure valid tip names at start of process

# outsider::module_install('dombennett/om..raxml')

# Functions ----
raxml_run <- function(wd, files_to_send, arglist) {
  otsdr <- outsider::.outsider_init(repo = "dombennett/om..raxml", 
                                    cmd = "raxmlHPC-PTHREADS-SSE3",
                                    wd = wd, files_to_send = files_to_send,
                                    arglist = arglist)
  outsider::.run(otsdr)
}

# Vars ----
supermatrix_dir <- file.path('primates', 'supermatrix')
raxml_dir <- file.path('primates', 'raxml')
if (!dir.exists(raxml_dir)) {
  dir.create(raxml_dir)
}
supermatrix_files <- list.files(path = supermatrix_dir, pattern = '.fasta')
partition_files <- list.files(path = supermatrix_dir, pattern = '.txt')


# Run ----
for (i in seq_along(supermatrix_files)) {
  nm <- sub(pattern = '\\.fasta', replacement = '', x = supermatrix_files[[i]])
  files_to_send <- c(file.path(supermatrix_dir, supermatrix_files[[i]]),
                     file.path(supermatrix_dir, partition_files[[i]]))
  arglist <- c('-f', 'a', '-m', 'GTRGAMMA', '-T', '1', '-#', '10', '-x',
               sample(1:10000000, 1), '-n', nm, '-s', supermatrix_files[[i]],
               '-q', partition_files[[i]], '-p', sample(1:10000000, 1))
  success <- raxml_run(wd = raxml_dir, files_to_send = files_to_send,
                       arglist = arglist)
  if (success) {
    file.remove(file.path(raxml_dir, partition_files[[i]]))
    file.remove(file.path(raxml_dir, partition_files[[i]]))
  }
}
