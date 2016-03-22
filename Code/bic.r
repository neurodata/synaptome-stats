
kbic <- function(x){
###
###  Input: a list of kmeans objects for k=1:K.
###  Output: a vector of BIC_k
###
###  I'm pretty sure there is a better 
###  way to code this ...

    K <- length(x)
    bic <- foreach(i = 1:K, .combine=c)%do%{
    kvec <- x[[i]]
    k <- nrow(kvec$centers)
    d <- ncol(kvec$centers)
    n <- length(kvec$cluster)

    sigHat <- 1/(length(kvec$cluster)-i)*kvec$tot.withinss
    lj <- foreach(j = 1:length(kvec$size),.combine=sum)%do%{
        -(kvec$size[j])/2*log(2*pi) -
            (kvec$size[j]*ncol(kvec$centers))/2 *
            log(sigHat) - (kvec$size[j] - i)/2 + 
            kvec$size[j] * log(kvec$size[j]) - 
            kvec$size[j] * log(length(kvec$cluster))
        }
    
    pj <- ((k-1) + d*k + 1)
    bic <- lj - pj/2 * log(length(kvec$cluster))
    bic
    }
    return(bic)
}
