#' Transform data in various ways for easy downstream manipulation.
#'
#' @param x A matrix, data.frame, or data.table with the data of
#' interest.
#' @param type A character list of the transformations desired. See
#' 'Details'
#' @return A list of the transformed data as data.tables.
#' @details If there are negatives present, using any of log transforms
#' will generate an error. Types of transformations incldue "Zscore" -- 
#' transforms to z-scores, "rank" -- passes to ranks, "1e3" -- scales 
#' the data between $[0,1e3]$, "log10" -- checks for negatives
#' and then removes 0's and takes the $\log_10$ of the data, "log" --
#' same as "log10" but for base $e$.  "slog10" -- checks for negatives
#' and then removes 0's and takes $\log_10$, then zscores. "slog" --
#' same as "slog10" but for base $e$. And "all" performs all of the
#' above.
#'
#' @examples
#' x <- runif(100)
#' y <- x + rnorm(100)
#' z <- x + max(y)
#' toyData <- cbind(x,y,z)
#' transData <- transforms(toyData,type=c("all"),
#' ties.method='average')
#'
transforms <- function(x, type = c('zscore'), ...) {

  type <- if (any(type == "all")) { 
    c("zscore", "rank", "1e3", "log10", "log", "slog10", "slog") 
  } else {
    type 
  }

  d <- as.data.table(x)
  out <- list(raw = d)
  
  if ("zscore" %in% type) {
    dZscore <- d[, lapply(.SD, scale, center = TRUE, scale = TRUE)]
    out[["Zscore"]] <- dZscore
  }
  
  if ("1e3" %in% type) {
    d01e3 <- d[, lapply(.SD, function(y) {
        (y - min(y))/(max(y) - min(y) * 1000)
    })]
    out[["d01e3"]] <- d01e3
  }
  
  if ("rank" %in% type) {
    dRank <- d[, lapply(.SD, rank, ...)]
    out[["Rank"]] <- dRank
  }
  
  negs <- any(d < 0)

  if ("log10" %in% type) {
    if (negs) stop("There are negatives!") 
    zs <- apply(d, 1, function(row) {
        any(row == 0)
    }) ## logical of which rows have 0's present.

    ## Select only rows greater than zero 
    dLog10 <- d[which(!zs == TRUE), lapply(.SD, log10)]
    out[["log10"]] <- dLog10
  }

  if ("log" %in% type) {
    if (negs) stop("There are negatives!") 
    zs <- apply(d, 1, function(row) {
        any(row == 0)
    }) ## logical of which rows have 0's present.

    ## Select only rows greater than zero 
    dLog <- d[which(!zs == TRUE), lapply(.SD, log)]
    out[["log"]] <- dLog
  }

  if ("slog10" %in% type) {
    if (negs) stop("There are negatives!") 
    zs <- apply(d, 1, function(row) {
        any(row == 0)
    }) ## logical of which rows have 0's present.

    ## Select only rows greater than zero 
    dlog10 <- d[which(!zs == TRUE), lapply(.SD, log10)]
    dSlog10 <- data.table(scale(dlog10, center=TRUE, scale=TRUE))

    out[["slog10"]] <- dSlog10
  }

  if ("slog" %in% type) {
    if (negs) stop("There are negatives!") 
    zs <- apply(d, 1, function(row) {
        any(row == 0)
    }) ## logical of which rows have 0's present.

    ## Select only rows greater than zero 
    dlog <- d[which(!zs == TRUE), lapply(.SD, log)]

    dSlog <- data.table(scale(dlog, center=TRUE, scale=TRUE))
    out[["slog"]] <- dSlog
  }
  
  return(out)
}

