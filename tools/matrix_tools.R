
seqs_read <- function(fl) {
  all_data <- readLines(fl)
  seqs <- list()
  for (i in seq_along(all_data)) {
    bit <- all_data[[i]]
    if (grepl(pattern = '^>', x = bit)) {
      nm <- sub(pattern = '^>', '', x = bit)
      seqs[[nm]] <- ''
    } else {
      bit <- strsplit(x = bit, split = '')[[1]]
      seqs[[nm]] <- c(seqs[[nm]], bit)
    }
  }
  class(seqs) <- 'sequences'
  seqs
}

print.sequences <- function(x) {
  max_ntips <- 10
  max_nclmns <- 50
  ntips <- ifelse(length(x) >= max_ntips, max_ntips, length(x))
  nclmns <- ifelse(length(x[[1]]) >= max_nclmns, max_nclmns, length(x[[1]]))
  part <- x[seq_len(ntips)]
  part <- lapply(X = part, FUN = function(x) paste0(x[seq_len(nclmns)],
                                                    collapse = ''))
  cat('Sequences: [', length(x), '] tips, [', length(x[[1]]), '] bps\n',
      sep = '')
  for (i in seq_along(part)) {
    cat('.... $', names(part)[i], ': ', part[[i]], ' ...\n', sep = '')
  }
  cat('.... ....\n')
}

gapmatrix_get <- function(seqs) {
  calc <- function(seq) {
    seq == '-'
  }
  nbps <- length(seqs[[1]])
  res <- t(vapply(X = seqs, FUN = calc, FUN.VALUE = logical(nbps)))
  class(res) <- 'gapmatrix'
  res
}

seqs_select <- function(seqs_list, nms) {
  seqs_get <- function(i) {
    seqs <- seqs_list[[i]]
    nbp <- length(seqs[[1]])
    filler <- rep('-', nbp)
    res <- vector(mode = 'list', length = length(nms))
    names(res) <- nms
    present <- names(seqs)[names(seqs) %in% nms]
    res[present] <- seqs[present]
    absents <- nms[!nms %in% names(seqs)]
    for (absent in absents) {
      res[[absent]] <- filler
    }
    class(res) <- 'sequences'
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
    seqs <- lapply(X = seqs, FUN = function(x) x[keep_clmns])
    class(seqs) <- 'sequences'
    seqs
  }
  seqs_list <- lapply(X = seqs_list, FUN = calc)
  nbps <- vapply(X = seqs_list, FUN = function(x) length(x[[1]]),
                 FUN.VALUE = integer(1))
  seqs_list[nbps >= min_nbps]
}

smatrix_get <- function(seqs_list) {
  stick_together <- function(i) {
    seqs <- lapply(X = seqs_list, FUN = '[[', i)
    seqs <- lapply(X = seqs, FUN = paste0, collapse = '')
    seq <- paste0(seqs, collapse = '')
  }
  nbps <- vapply(X = seqs_list, FUN = function(x) length(x[[1]]),
                 FUN.VALUE = integer(1))
  res <- lapply(seq_along(seqs_list[[1]]), stick_together)
  names(res) <- names(seqs_list[[1]])
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


