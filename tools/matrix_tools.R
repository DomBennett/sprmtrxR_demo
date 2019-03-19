
seqs_read <- function(fl) {
  all_data <- readLines(fl)
  seqs <- list()
  for (i in seq_along(all_data)) {
    bit <- all_data[[i]]
    if (grepl(pattern = '^>', x = bit)) {
      nm <- sub(pattern = '^>', '', x = bit)
      seqs[[nm]] <- NULL
    } else {
      bit <- strsplit(x = bit, split = '')[[1]]
      seqs[[nm]] <- c(seqs[[nm]], bit)
    }
  }
  lngths <- vapply(X = seqs, FUN = length, FUN.VALUE = integer(1))
  if (all(lngths != lngths[[1]])) {
    stop('Not an alignment: sequences have different lengths.')
  }
  res <- matrix(unlist(seqs), ncol = length(seqs[[1]]), nrow = length(seqs),
                byrow = TRUE)
  rownames(res) <- names(seqs)
  class(res) <- 'sequence_alignment'
  res
}

print.sequence_alignment <- function(x) {
  max_ntips <- 10
  max_nclmns <- 50
  ntips <- ifelse(nrow(x) >= max_ntips, max_ntips, length(x))
  nclmns <- ifelse(ncol(x) >= max_nclmns, max_nclmns, ncol(x))
  part <- x[seq_len(ntips), seq_len(nclmns)]
  cat('Sequence alignment: [', nrow(x), '] tips, [', ncol(x),
      '] bps\n',
      sep = '')
  for (i in seq_len(ntips)) {
    cat('.... $', rownames(part)[i], ': ', part[i, ], ' ...\n', sep = '')
  }
  cat('.... ....\n')
}

gapmatrix_get <- function(seqs) {
  res <- seqs == '-'
  rownames(res) <- rownames(seqs)
  class(res) <- 'gapmatrix'
  res
}

seqs_select <- function(seqs_list, nms) {
  seqs_get <- function(i) {
    seqs <- seqs_list[[i]]
    nbp <- ncol(seqs)
    res <- matrix(data = '-', nrow = length(nms), ncol = nbp)
    rownames(res) <- nms
    present <- rownames(seqs)[rownames(seqs) %in% nms]
    res[present, ] <- seqs[present, ]
    class(res) <- 'sequence_alignment'
    res
  }
  res <- lapply(seq_along(seqs_list), seqs_get)
  names(res) <- names(seqs_list)
  res
}

seqs_filter <- function(seqs_list, cutoff = 0.9, min_nbps = 200) {
  calc <- function(seqs) {
    gapmatrix <- gapmatrix_get(seqs = seqs)
    pclmn <- 1 - (colSums(gapmatrix)/nrow(gapmatrix))
    keep_clmns <- pclmn >= cutoff
    seqs <- seqs[ ,keep_clmns]
    class(seqs) <- 'sequence_alignment'
    seqs
  }
  seqs_list <- lapply(X = seqs_list, FUN = calc)
  nbps <- vapply(X = seqs_list, FUN = ncol, FUN.VALUE = integer(1))
  seqs_list[nbps >= min_nbps]
}

smatrix_get <- function(seqs_list) {
  stick_together <- function(i) {
    res <- unlist(lapply(X = seqs_list, FUN = function(x) x[i, ]))
    paste0(res, collapse = '')
  }
  nbps <- vapply(X = seqs_list, FUN = ncol, FUN.VALUE = integer(1))
  ntips <- nrow(seqs_list[[1]])
  res <- lapply(seq_len(ntips), stick_together)
  names(res) <- rownames(seqs_list[[1]])
  attr(res, 'genes') <- names(seqs_list)
  attr(res, 'nbps') <- nbps
  class(res) <- 'smatrix'
  res
}

print.smatrix <- function(x) {
  nmssng <- sum(vapply(X = gregexpr(pattern = '-', text = x), FUN = length,
                       FUN.VALUE = integer(1)))
  total_nbps <- sum(attr(x, 'nbps')) * length(x)
  pmssng <- signif(x = nmssng * 100/total_nbps, digits = 2)
  cat('Smatrix:\n', '.... [', length(attr(x, 'genes')), '] genes\n',
      '.... [', length(x), '] tips\n',
      '.... [', sum(attr(x, 'nbps')), '] bps long\n',
      '.... [', pmssng, '%] gaps\n', sep = '')
}

drop_tips <- function(smatrix, cutoff = 0.5) {
  nmissing <- vapply(X = gregexpr(pattern = '-', text = smatrix),
                     FUN = length, FUN.VALUE = integer(1))
  nbps <- sum(attr(smatrix, 'nbps'))
  pmissing <- 1 - (nmissing/nbps)
  keep_tips <- pmissing >= cutoff
  res <- smatrix[keep_tips]
  attr(res, 'genes') <- attr(smatrix, 'genes')
  attr(res, 'nbps') <- attr(smatrix, 'nbps')
  class(res) <- 'smatrix'
  res
}

smatrices_get <- function(groups, alignments, column_cutoff = .5,
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
    if (length(seqs_list) == 0) {
      next
    }
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
  class(res) <- 'smatrices'
  res
}

print.smatrices <- function(x) {
  cat('Smatrices object of [', length(x), ']', sep = '')
}
