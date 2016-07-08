\name{bhkmpp}
\alias{bhkmpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Binary-Hierarchical K-means++
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Runs a hierarchical K-means++ for $k=2$ on $n$ levels.
}
\usage{
bhkmpp(x, blevels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
The data.
}
  \item{blevels}{
%%     ~~Describe \code{blevels} here~~
The number of levels.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A data.table of cluster labels is returned where column $i$ holds
the labels for cluster $i$.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Jesse Leigh Patsolic
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, blevels) 
{
    require(data.table)
    k0 <- kmpp(x, 2)
    L <- data.table(lv1 = k0$cluster)
    for (y in 2:blevels) {
        L[[y]] <- 0L
    }
    colnames(L) <- paste0("lv", 1:blevels)
    for (j in 1:(blevels - 1)) {
        for (i in sort(unique(L[[j]]))) {
            set.seed(13)
            kv <- kmpp(x[L[[j]] == i, ], k = 2)
            if (i != 1) {
                kv$cluster <- as.integer(kv$cluster + 2 * (i - 
                  1))
            }
            L[L[[j]] == i, ][[j + 1]] <- kv$cluster
        }
    }
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
