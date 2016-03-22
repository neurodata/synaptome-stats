clusterFraction <- function(kvec){
## Given a K-means object return a labeling that includes
##  the cluster number and fraction of elements within 
##  that cluster
##
## Args:
##   kvec: a k-means object or vector of cluster labels
##
## Return:
##   cs: a vector of length ncol(kvec$centers) formatted
##          e.g. "C1 0.45".  That is Cluster 1, 45% of total.

require(foreach)

if(class(kvec) == "kmeans"){
    labs <- kvec$cluster
    } else {
        labs <- kvec
    }

cs <- paste(paste0("C",1:max(labs)," "),
   foreach(i =1:max(labs),.combine=c)%do%{(sprintf("%.2g",sum(labs==i)/length(labs)))})

return(cs)
}
