#' Feature Exploration: (Column Centric)
#'
#' @param x A single data.table or a list of data.tables of features.
#' @param t98 A boolean which determines if the data are truncated
#' to the inner 98%. 
#' @return A list of various exploratory plots.
#'
#' @examples
#' x <- runif(100)
#' y <- x + rnorm(100)
#' z <- x + max(y)
#' toyData <- cbind(x,y,z)
#'
featureEx <- 
function(x, t98=NULL) {
  d <- as.data.table(x)

  if (!is.null(t98)) {
    ### Transform data here 
  }

    
}
