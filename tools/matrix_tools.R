overhangs_drop <- function(als, cutoff = .75) {
  for (i in seq_along(als)) {
    gps_mtrx <- matrix(FALSE, ncol = nchar(als[[i]][[1]]),
                       nrow = length(als[[i]]))
    for(j in seq_along(als[[i]])) {
      pull <- gregexpr('-', als[[i]][[j]])[[1]]
      if(pull[1] == -1) {
        next
      }
      gps_mtrx[j, pull] <- TRUE
    }
    prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
    ovrlppng <- which(prp_mssng > cutoff)
    strt <- ovrlppng[1]
    end <- ovrlppng[length(ovrlppng)]
    for(j in seq_along(als[[i]])) {
      als[[i]][[j]] <- substr(als[[i]][[j]], start = strt, stop = end)
    }
  }
  als
}

split_up <- function(al, cutoff = .75, min_length = 500) {
  gps_mtrx <- matrix(FALSE, ncol = nchar(al[[1]]), nrow = length(al))
  for(i in seq_along(al)) {
    pull <- gregexpr('-', al[[i]])[[1]]
    gps_mtrx[i, pull] <- TRUE
  }
  prp_mssng <- 1 - (colSums(gps_mtrx)/nrow(gps_mtrx))
  strts <- which((prp_mssng[-1] > cutoff) &
                   (prp_mssng[-1*length(prp_mssng)] < cutoff)) + 1
  ends <- which((prp_mssng[-1] < cutoff) &
                  (prp_mssng[-1*length(prp_mssng)] > cutoff))
  lngths <- ends - strts
  pull <- lngths > min_length
  strts <- strts[pull]
  ends <- ends[pull]
  new_als <- vector('list', length = length(strts))
  for(i in seq_along(strts)) {
    strt <- strts[i]
    end <- ends[i]
    new_al <- vector('list', length = length(al))
    for(j in seq_along(al)) {
      new_al[[j]] <- substr(al[[j]], start = strt, stop = end)
    }
    names(new_al) <- names(al)
    new_als[[i]] <- new_al
  }
  new_als
}

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

seqs_write <- function(sqs, fl) {
  fasta <- ''
  for(i in seq_along(sqs)) {
    sq <- sqs[[i]]
    dfln <- names(sqs)[[i]]
    fasta <- paste0(fasta, '>', dfln, '\n', sq, '\n\n')
  }
  cat(fasta, file = fl)
}

partition <- function(lngths, fl) {
  gene <- strt <- 1
  prttn_txt <- ''
  for(lngth in lngths) {
    end <- lngth + strt - 1
    prttn_txt <- paste0(prttn_txt, 'DNA, gene', gene, ' = ', strt, '-', end,
                        '\n')
    strt <- end + 1
    gene <- gene + 1
  }
  cat(prttn_txt, file = fl)
}

smatrix_get <- function(nms, alignments, nbps) {
  seqs_get <- function(i) {
    nbp <- nbps[[i]]
    algnmnt <- alignments[[i]]
    filler <- rep('-', nbp)
    res <- vector(mode = 'list', length = length(nms))
    names(res) <- nms
    present <- names(algnmnt)[names(algnmnt) %in% nms]
    res[present] <- algnmnt[present]
    absents <- nms[!nms %in% names(algnmnt)]
    for (absent in absents) {
      res[[absent]] <- filler
    }
    class(res) <- 'sequences'
    res
  }
  seqs <- lapply(seq_along(alignments), seqs_get)
  seqs <- lapply(seqs, paste0, collapse = '')
  res <- list('seqs' = seqs, 'genes' = names(alignments), 'nbases' = nbases,
              tips = nms)
  class(res) <- 'smatrix'
  res
}

print.smatrix <- function(x) {
  nmssng <- sum(vapply(X = gregexpr(pattern = '-', text = x$seqs), FUN = length,
                       FUN.VALUE = integer(1)))
  '-' == x$seqs[[1]]
  
  pmssng <- signif(x = nmssng * 100/(sum(x$nbases) * length(x$seqs)),
                   digits = 2)
  cat('Smatrix:\n', '.... [', length(x$genes), '] genes\n',
      '.... [', length(x$tips), '] tips\n',
      '.... [', sum(x$nbases), '] bps long\n',
      '.... [', pmssng, '%] gaps\n', sep = '')
}

gapmatrix_get <- function(smatrix) {
  calc <- function(seq) {
    seq == '-'
  }
  nbps <- length(smatrix[[1]])
  res <- t(vapply(X = smatrix, FUN = calc, FUN.VALUE = logical(nbps)))
  class(res) <- 'gapmatrix'
  res
}

drop_columns <- function(smatrix, cutoff = 0.9) {
  calc <- function(seq) {
    seq[keep_clmns]
  }
  gapmatrix <- gapmatrix_get(smatrix = smatrix)
  pclmn <- 1 - (colSums(gapmatrix)/ncol(gapmatrix))
  keep_clmns <- pclmn >= cutoff
  smatrix <- lapply(X = smatrix, FUN = calc)
  smatrix
}

drop_tips <- function(smatrix, cutoff = 0.5) {
  gapmatrix <- gapmatrix_get(smatrix = smatrix)
  ptips <- 1 - (rowSums(gapmatrix)/ncol(gapmatrix))
  keep_tips <- ptips >= cutoff
  smatrix[keep_tips]
}
