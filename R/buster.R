#'
#' Buster
#'
#' Performs bagging on the hierarchical clusters
#'
#' @param dist A distance object
#' @param n The number of times the data should be resampled
#' @param k The number of clusters
#' @param size The percentage of the original data in each resample
#' @param method The linkage method to be passed to hclust
#' @param outlier.th The threshold below which observations should be discarded
#' @return An an object of class buster which includes an hclust object on the co-ocurrences, an hclust object on the original distance measure and an evalution of the stability of the observations 
#' for model objects specific to the type of classifier. 
#' @author Simon Raper
#' @examples
#' 
#' us.dist<-dist(USArrests)
#' bhc<-buster(us.dist, n=250, k=5, size=0.66, method='ward', outlier.th=0.1)
#' plot(bhc)
#' 
#' #Identifies the states in the California cluster as being very voltile. We expect this cluster to break up
#' 
#' #Simple test
#' 
#' #First look at how it picks out the borderline cases
#' 
#' x.1<-rnorm(50, 10, 3)
#' y.1<-rnorm(50, 10, 3)
#' x.2<-rnorm(50, 12, 3)
#' y.2<-rnorm(50, 15, 4)
#' x.3<-rnorm(50, 10, 3)
#' y.3<-rnorm(50, 13, 4)
#' 
#' test.data<-data.frame(group=rep(1:3, each=50), x=c(x.1, x.2, x.3), y=c(y.1, y.2, y.3))
#' names<-c(paste0("group 1: ", 1:50), paste0("group 2: ", 1:50), paste0("group 3: ", 1:50))
#' rownames(test.data)<-names
#' dist<-dist(test.data[,-1])
#' 
#' bhc<-buster(dist, n=200, k=3, size=0.66, method='ward', outlier.th=0.9)
#' 
#' plot(bhc)
#' 
#' #Append the assigned clusters
#' graph.data<-cbind(test.data, bhc$obs.eval[order(bhc$obs.eval$obs.ind, decreasing=FALSE),])
#' plot(graph.data$x, graph.data$y, xlim=c(0,30), ylim=c(0, 30), pch = graph.data$group, col=graph.data$cluster+1)
#' 
#' max.co<-bhc$obs.eval$max.co[order(bhc$obs.eval$obs.ind)]
#' alpha<-(max.co-min(max.co))/(max(max.co)-min(max.co))
#' cols <- hsv(0,0,0,alpha)
#' plot(graph.data$x, graph.data$y, xlim=c(0,30), ylim=c(0, 30), pch = 19, col=cols)



buster<-function(dist, n=100, k, size=0.66, method='ward', outlier.th = 0.7) {
  
  #Constants
  dist.m<-as.matrix(dist)
  nd<-nrow(dist.m)
  size<-round(size*nd)
  
  #Store the results of n clusterings
  clus.bs<-NULL
  for (i in 1:n) {
    bs.ind<-sample.int(nd, size, replace=FALSE)
    samp.dist<-as.dist(dist.m[bs.ind, bs.ind])
    hc <- hclust(samp.dist, method)
    ct<-cutree(hc, k)
    add<-data.frame(iter=rep(i, size), ind=names(ct), cluster=ct)
    clus.bs<-rbind(clus.bs, add)
  }
    
  #Work out co-occurence in clusters for the observations
  m<-merge(clus.bs, clus.bs, by=c("iter", "cluster"))
  
  d<-dcast(m, ind.x~ind.y, length)
  rownames(d)<-d[,1]
  d<-as.matrix(d[rownames(dist.m),rownames(dist.m)])
  
  #Work out co-appearances for the observations
  m2<-merge(clus.bs, clus.bs, by="iter")
  
  d2<-dcast(m2, ind.x~ind.y, length)
  rownames(d2)<-d2[,1]
  d2<-as.matrix(d2[rownames(dist.m),rownames(dist.m)])

  #Create the distance matrix
  dm<-d/d2
  disim<-as.dist(1-dm)
  
  browser()
  
  #Work out which observations are promiscuous
  dm.rep<-dm
  diag(dm.rep)<-0
  most<-apply(dm.rep,1,var)
  #exclude<-which(most<outlier.th)
  include<-which(most>=outlier.th)
  
  #Dendrogram with promiscuous observations labelled
  
  h<-hclust(disim, "ward")
  
  #Remove them from the distance matrix
  
  dm<-d/diag(d)
  disim<-as.dist(1-dm)
  
  #Run hierarchical clustering
  h.co<-hclust(disim, method)
  
  eval<-data.frame(obs.ind=1:nd, obs.names=names(most), max.co = most, exclude=most<=outlier.th)
  clus<-data.frame(obs.names=h$labels, cluster = cutree(h, k))
  m<-merge(clus, eval, by="obs.names", all.x=TRUE)
  m$cluster[m$exclude==TRUE]<-0

  buster<-list(bhclust=h.co, hclust=h, obs.eval=m)
  
  class(buster)<-"buster"
  
  return(buster)
  
}
