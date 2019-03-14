nms_from_fasta <- function(flpth) {
  lines <- readLines(con = flpth)
  nms <- lines[grepl(pattern = '^>', x = lines)]
  sub(pattern = '^>', replacement = '', nms)
}

alignment_info_get <- function(flpth) {
  nms <- nms_from_fasta(flpth = flpth)
  res <- list(nms = nms, ntips = length(nms), flpth = flpth)
  class(res) <- 'alignment_info'
  res
}

nms_from_alignments <- function(flpths) {
  calc <- function(flpth) {
    algnmnt <- alignment_info_get(flpth = flpth)
    algnmnt[['nms']]
  }
  nms <- lapply(X = flpths, FUN = calc)
  unique(unlist(nms))
}

print.alignment_info <- function(x) {
  cat('Alignment (', outsider:::stat(x$ntips), ') tips.\n', sep = '')
}

nm_match <- function(algnmnt_nms, tree_nms, algnmnt_pttrns = algnmnt_nms,
                     tree_pttrns = tree_nms, max_dist = .1) {
  # Checks
  if (any(duplicated(algnmnt_nms)) | any(duplicated(tree_nms))) {
    stop('Duplicated names.')
  }
  if (length(algnmnt_nms) != length(algnmnt_pttrns) |
      length(tree_nms) != length(tree_pttrns)) {
    stop('`*_nms` and `*_pttrns` must be same lengths.')
  }
  calc <- function(i) {
    # Calc prop. Levenshtein distance, select tree name with lowest distance
    algnmnt_pttrn <- algnmnt_pttrns[[i]]
    dists <- adist(x = algnmnt_pttrn, y = tree_pttrns, partial = TRUE)[1, ]
    pdists <- dists/nchar(algnmnt_pttrn)
    pssbls <- which(pdists < max_dist)
    npssbls <- length(pssbls)
    if (npssbls > 1) {
      res <- tree_nms[pssbls[which.min(pdists[pssbls])]]
    } else if (npssbls == 1) {
      res <- tree_nms[pssbls]
    } else {
      res <- ''
    }
    res
  }
  matched_treenms <- vapply(X = seq_along(algnmnt_nms), FUN = calc,
                            FUN.VALUE = character(1))
  unmatched <- algnmnt_nms[matched_treenms == '']
  matched <- algnmnt_nms[matched_treenms != '']
  matched_treenms <- matched_treenms[matched_treenms != '']
  res <- list('alignment' = matched, 'tree' = matched_treenms,
              'unmatched' = unmatched)
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
