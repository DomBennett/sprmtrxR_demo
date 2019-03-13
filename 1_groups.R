# Functions ----
nms_from_fasta <- function(flpth) {
  lines <- readLines(con = flpth)
  nms <- lines[grepl(pattern = '^>', x = lines)]
  sub(pattern = '^>', replacement = '', nms)
}

alignment_get <- function(flpth) {
  nms <- nms_from_fasta(flpth = flpth)
  res <- list(nms = nms, ntips = length(nms), flpth = flpth)
  class(res) <- 'alignment'
  res
}

nms_from_alignments <- function(flpths) {
  calc <- function(flpth) {
    algnmnt <- alignment_get(flpth = flpth)
    algnmnt[['nms']]
  }
  nms <- lapply(X = flpths, FUN = calc)
  unique(unlist(nms))
}

print.alignment <- function(x) {
  cat('Alignment (', outsider:::stat(x$ntips), ') tips.\n', sep = '')
}

nm_match <- function(algnmnt_nms, tree_nms, max_dist = .1) {
  if (any(duplicated(algnmnt_nms)) | any(duplicated(tree_nms))) {
    stop('Duplicated names.')
  }
  calc <- function(algnmnt_nm) {
    # Calc prop. Levenshtein distance, select tree name with lowest distance
    dists <- adist(x = algnmnt_nm, y = tree_nms, partial = TRUE)[1, ]
    pdists <- dists/nchar(algnmnt_nm)
    pssbls <- which(pdists < max_dist)
    npssbls <- length(pssbls)
    if (npssbls > 2) {
      res <- tree_nms[pssbls[which.min(pdists[pssbls])]]
    } else if (npssbls == 1) {
      res <- tree_nms[pssbls]
    } else {
      res <- ''
    }
    res
  }
  res <- vapply(X = algnmnt_nms, FUN = calc, FUN.VALUE = character(1))
  nms <- res[res != '']
  matched <- names(nms)
  matched <- unname(matched)
  unmatched <- names(res[res == ''])
  algnmnts_nms_wdups <- matched[duplicated(nms)]
  if (length(algnmnts_nms_wdups) > 0) {
    msg <- paste0('Alignment names are matching to the same tree name(s):\n',
                  paste0(algnmnts_nms_wdups, collapse = ', '))
    warning(msg)
  }
  res <- list('alignment' = matched, 'tree' = nms, 'unmatched' = unmatched)
  class(res) <- 'matched_names'
  res
}

print.matched_names <- function(x) {
  format_text <- function(y) {
    n <- ifelse(length(y) > 3, 3, length(y))
    paste0(paste0(y[1:n], collapse = ', '), ' ...')
  }
  cat('[', length(x$alignment), '] names matched:\n',
      '... from `$alignment` ', format_text(x$alignment),
      '\n... to `$tree` ', format_text(x$tree),
      '\n[', length(x$unmatched), '] unmatched.\n', sep = '')
}

groups_get <- function(matched_nms, tree, max_size, min_size) {
  group_get <- function(tree, groups) {
    ptids <- treeman::getNdSlt(tree = tree, slt_nm = 'ptid', id = tree@root)
    ptids <- ptids[!ptids %in% tree@tips]
    for (ptid in ptids) {
      subtree <- treeman::getSubtree(tree = tree, id = ptid)
      pssbls <- unname(nms[nms %in% subtree@tips])
      size <- length(pssbls)
      if (size <= max_size & size >= min_size) {
        groups[[ptid]] <- pssbls
      } else {
        if (subtree@ntips > min_size) {
          groups <- group_get(tree = subtree, groups = groups)
        }
      }
    }
    groups
  }
  nms <- matched_nms$tree
  alignment_nms <- matched_nms$alignment
  groups <- group_get(tree = tree, groups = list())
  res <- lapply(X = groups, function(x) alignment_nms[nms %in% x])
  unmatched <- c(matched_nms[['unmatched']],
                 alignment_nms[!alignment_nms %in% unlist(res)])
  res[['unmatched']] <- unmatched
  class(res) <- "groups"
  res
}

print.groups <- function(x) {
  nmatched <- vapply(X = x, FUN = length, integer(1))
  nunmatched <- nmatched[['unmatched']]
  nmatched <- sum(nmatched[names(nmatched) != 'unmatched'])
  cat('Tips group object of [', length(x), '] elements:\n', '... [', nmatched,
      '] tips in mono\n... [', nunmatched, '] tips in super',
      sep = '')
}

# Vars ----
max_size <- 10
min_size <- 5
input_dir <- file.path('primates', 'alignments')
alfls <- file.path(input_dir, list.files(path = input_dir, pattern = '.fasta'))

# Data ----
tree <- treeman::readTree(file = file.path('0_data', 'primate.tre'))
tree_nms <- tree@tips
algnmnt_nms <- nms_from_alignments(flpths = alfls)
# TODO add patterns
matched_nms <- nm_match(algnmnt_nms = algnmnt_nms, tree_nms = tree_nms,
                        max_dist = .15)
groups <- groups_get(matched_nms = matched_nms, tree = tree,
                     max_size = max_size, min_size = min_size)
groups





