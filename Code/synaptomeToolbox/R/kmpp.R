#' K-means++ initialization followed by K-means.
#'
#' @param x A matrix, data.frame, or data.table 
#' @param k An integer denoting how many centers.
#' @return A kmeans object. 
#' @examples
#' data(iris)
#' kmpp(iris[,1:4],k=3)
kmpp <- function(x, k = 2) {
    
    la2 <- function(a, b) {
        m1 <- matrix(as.numeric(a), nrow = dim(b)[1], ncol = dim(b)[2], byrow = TRUE)
        
        if (all(dim(m1) == dim(b))) {
            dif <- apply((m1 - b)^2, 1, sum)
            return(dif)
        } else {
            stop("Wrong dimensions!")
        }
    }
    
    ### Step 1a
    s1 <- cp <- sample(nrow(x), 1)
    
    ind <- c(cp)
    
    centers <- list()
    
    centers[[1]] <- x[s1, ]
    
    p <- la2(x[s1, ], x)
    P <- p/sum(p)  ## P = D(x')^2/( ∑_x\inX  D(x)^2)
    
    ### Step 1b-c
    for (i in 2:k) {
        # set.seed(2^13)
        cp <- sample(nrow(x), 1, prob = P)  ## c'
        centers[[i]] <- x[cp, ]
        ind <- c(ind, cp)
        
        p <- apply(cbind(la2(x[cp, ], x), p), 1, min)  #apply(cbind(p,ds[,cp]),1,min)
        P <- p/sum(p)
    }
    
    Ctrs <- Reduce("rbind", centers)
    
    kObj <- kmeans(x, centers = Ctrs)
    
    return(kObj)
}
