#' K-Means++ initialization to K-means
#'
#' INPUTS:
#' x: The data formatted as a matrix.
#' k: The number of clusters desired. \code{k < nrows(x)}
#' 
#'
#' REFERENCE:
#' Arthur, David and Vassilvitskii, Sergei k-means++: The Advantages of Careful Seeding
#'


kmpp <- function(x,k=2){

la2 <- function(a,b){
        m1 <- matrix(as.numeric(a),nrow=dim(b)[1],ncol=dim(b)[2], byrow=TRUE)
    
        if(all(dim(m1) == dim(b))){
            dif <- apply((m1-b)^2,1,sum)    
        return(dif)
        } else {
            stop("Wrong dimensions!")}
        }

### Step 1a
#set.seed(2^13)
s1 <- cp <- sample(nrow(x),1)

ind <- c(cp)

centers <- list()

centers[[1]] <- x[s1,]

p <- la2(x[s1,],x)
P <- p/sum(p) ## P = D(x')^2/( ∑_x\inX  D(x)^2)

### Step 1b-c
for( i in 2:k ) {
    #set.seed(2^13)
    cp <- sample(nrow(x),1,prob=P) ## c'
    centers[[i]] <-x[cp,]
    ind <- c(ind,cp)
    
    p <- apply(cbind(la2(x[cp,],x),p),1,min) #apply(cbind(p,ds[,cp]),1,min)
    P <- p/sum(p)
    }

Ctrs <- Reduce('rbind',centers)

kObj <- kmeans(x, centers=Ctrs)

return(kObj)
} ### END FUNCTION


bhkmpp <- function(x,blevels){
    require(data.table)

    k0 <- kmpp(x,2)
    L <- data.table(lv1=k0$cluster)

    for(y in 2:blevels){ L[[y]] <- 0L }
    colnames(L) <- paste0("lv", 1:blevels)
    
    for(j in 1:(blevels-1)){
    for(i in sort(unique(L[[j]]))){
        ## Set seed to keep cluster labels somewhat consistant 
        set.seed(13)
        kv <- kmpp(x[L[[j]] == i,], k = 2)
        if(i != 1){kv$cluster <- as.integer(kv$cluster + 2*(i-1))}
        L[L[[j]]==i,][[j+1]] <- kv$cluster
        }
    }
    return(L)
}




