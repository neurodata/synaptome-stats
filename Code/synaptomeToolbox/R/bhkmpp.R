bhkmpp <- function(x, blevels) {
    
    k0 <- kmpp(x, 2)
    L <- data.table(lv1 = k0$cluster)
    
    for (y in 2:blevels) {
        L[[y]] <- 0L
    }
    colnames(L) <- paste0("lv", 1:blevels)
    
    for (j in 1:(blevels - 1)) {
        for (i in sort(unique(L[[j]]))) {
            ## Set seed to keep cluster labels somewhat consistant
            set.seed(13)
            kv <- kmpp(x[L[[j]] == i, ], k = 2)
            if (i != 1) {
                kv$cluster <- as.integer(kv$cluster + 2 * (i - 1))
            }
            L[L[[j]] == i, ][[j + 1]] <- kv$cluster
        }
    }
    return(L)
}
