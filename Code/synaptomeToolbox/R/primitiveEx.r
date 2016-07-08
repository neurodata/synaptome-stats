#' Primitive Exploration 
#'
#' @param x A list of data.frames or data.tables of primitives
#' @return A list of various plots using dimension reduction. 
#' @details Details to come. 
#' @examples
#' x <- runif(100)
#' y <- x + rnorm(100)
#' z <- x + max(y)
#' toyData <- cbind(x,y,z)
#' ### Finish examples
primitiveEx <- function(x, coi) {
  x <- dF

  cors <- lapply(x, cor)

  par(mfrow = c(length(cors),1))
  for (i in cors) {
    corrplot(i, method='color', tl.col='black', tl.cex=0.75)
  }

}
