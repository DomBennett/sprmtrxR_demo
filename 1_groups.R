# Functions ----
tipnms_from_fasta <- function(flpth) {
  lines <- readLines(con = flpth)
  tipnms <- lines[grepl(pattern = '^>', x = lines)]
  sub(pattern = '^>', replacement = '', tipnms)
}

alignment_get <- function(flpth) {
  tipnms <- tipnms_from_fasta(flpth = flpth)
  res <- list(tipnms = tipnms, ntips = length(tipnms), flpth = flpth)
  class(res) <- 'alignment'
  res
}

tipnms_from_alignments <- function(flpths) {
  calc <- function(flpth) {
    algnmnt <- alignment_get(flpth = flpth)
    algnmnt[['tipnms']]
  }
  tipnms <- lapply(X = flpths, FUN = calc)
  unique(unlist(tipnms))
}

print.alignment <- function(x) {
  cat('Alignment (', outsider:::stat(x$ntips), ') tips.\n', sep = '')
}

tipnms <- function(x) {
  UseMethod('tipnms', x)
}

# Vars ----
input_dir <- file.path('primates', 'alignments')
alfls <- file.path(input_dir, list.files(path = input_dir, pattern = '.fasta'))

# Data ----
tree <- ape::read.tree(file = file.path('0_data', 'primate.tre'))
tree_tipnms <- tree$tip.label
tree_tipnms <- sub(pattern = '_', replacement = ' ', x = tree_tipnms)
algnmnt_tipnms <- tipnms_from_alignments(flpths = alfls)
matched_tipnms <- tipnm_match(algnmnt_tipnms=algnmnt_tipnms,
                              tree_tipnms=tree_tipnms)



tipnm_match <- function(algnmnt_tipnms, tree_tipnms, max_dist = .1) {
  if (any(duplicated(algnmnt_tipnms)) | any(duplicated(tree_tipnms))) {
    stop('Duplicated names.')
  }
  calc <- function(algnmnt_tipnm) {
    # Calc prop. Levenshtein distance, select tree name with lowest distance
    dists <- adist(x = algnmnt_tipnm, y = tree_tipnms, partial = TRUE)[1, ]
    pdists <- dists/nchar(algnmnt_tipnm)
    pssbls <- which(pdists < max_dist)
    npssbls <- length(pssbls)
    if (npssbls > 2) {
      res <- tree_tipnms[pssbls[which.min(pdists[pssbls])]]
    } else if (npssbls == 1) {
      res <- tree_tipnms[pssbls]
    } else {
      res <- ''
    }
    res
  }
  res <- vapply(X = algnmnt_tipnms, FUN = calc, FUN.VALUE = character(1))
  res_womssng <- res[res != '']
  algnmnts_tipnms_wdups <- names(res_womssng)[duplicated(res_womssng)]
  if (length(algnmnts_tipnms_wdups) > 0) {
    msg <- paste0('Alignment names are matching to the same tree name(s):\n',
                  paste0(algnmnts_tipnms_wdups, collapse = ', '))
    warning(msg)
  }
  res
}
