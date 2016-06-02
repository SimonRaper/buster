#'
#' Buster
#'
#' Performs bagging on hierarchical clusters
#'
#' @param dist A distance object
#' @param n The number of times the data should be resampled
#' @param k The number of clusters
#' @param size The percentage of the original data in each resample
#' @param method The linkage method to be passed to hclust
#' @param pct.exc Setting this to x excludes or highlights the top x % observations ranked by instability. For example if we set this to 0.1 then the top decile of unstable observations are excluded or highlighted.
#' @param low.mem Setting this to true uses a slower but less memory intensive way of calculating the co-occurrences
#' @return An an object of class buster which includes
#' \itemize{
#'  \item{bhclust }{An hclust object on the co-ocurrences}
#'  \item{bhclust.ex }{An hclust object on the co-ocurrences exclduing unstable observations}
#'  \item{hclust }{An hclust object on the original distance measure}
#'  \item{obs.eval }{An evalution of the stability of the observations}
#' } 
#' @author Simon Raper
#' @examples
#' 
#' #Testing on the iris data set
#' iris.dist<-dist(iris[,1:4])
#' bhc<-buster(iris.dist, n=250, k=3, size=0.66, method='ward.D', pct.exc=0.07)
#' plot(bhc)
#' 
#' #We see the unstable observations in pink.
#' 
#' cluster<-bhc$obs.eval$cluster[order(bhc$obs.eval$obs.ind)]
#' plot(iris[,1:4], col=6-cluster, pch = rep(15:17, each=50))
#' 
#' #Another simple test
#' 
#' x.1<-rnorm(50, 10, 3)
#' y.1<-rnorm(50, 10, 3)
#' x.2<-rnorm(50, 20, 3)
#' y.2<-rnorm(50, 10, 3)
#' x.3<-rnorm(50, 13, 3)
#' y.3<-rnorm(50, 20, 3)
#' 
#' test.data<-data.frame(group=rep(1:3, each=50), x=c(x.1, x.2, x.3), y=c(y.1, y.2, y.3))
#' names<-c(paste0("group 1: ", 1:50), paste0("group 2: ", 1:50), paste0("group 3: ", 1:50))
#' rownames(test.data)<-names
#' dist<-dist(test.data[,-1])
#' 
#' bhc<-buster(dist, n=200, k=3, size=0.66, method='ward.D', pct.exc=0.1)
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

buster<-function(dist, n=100, k, size=0.66, method='ward.D', pct.exc=0.1) {
  
  #Constants
  dist.m<-as.matrix(dist)
  nd<-nrow(dist.m)
  size<-round(size*nd)
  
  #Set up matrices
  
  cs<-matrix(rep(0, nd*nd), nrow=nd)
  rownames(cs)<-rownames(dist.m)
  colnames(cs)<-colnames(dist.m)
  co<-matrix(rep(0, nd*nd), nrow=nd)
  rownames(co)<-rownames(dist.m)
  colnames(co)<-colnames(dist.m)
  
  #Store the results of n clusterings
  #clus.bs<-NULL
  for (i in 1:n) {
    bs.ind<-sample.int(nd, size, replace=FALSE) #Resample
    sampled<-as.numeric(1:nd %in% bs.ind) #Create a vector indicating whether or not a data point is included in the sample
    cs.inc<-sampled%*%t(sampled) #Increment for co-sampled
    cs<-cs+cs.inc
    samp.dist<-as.dist(dist.m[bs.ind, bs.ind])
    hc <- hclust(samp.dist, method)
    ct<-cutree(hc, k)
    ct<-data.frame(names=names(ct), ct=ct)
    names.list<-data.frame(names=rownames(dist.m))
    clusters<-merge(names.list, ct, by="names", all.x=TRUE)
    clusters$ct[is.na(clusters$ct)]<-0
    clusters$ct<-as.factor(clusters$ct)
    
    dummied<-model.matrix(~ct, data=clusters)
    rownames(dummied)<-clusters$names
    
    dummied[rowSums(dummied[,-1])>0,1]<-0
    dummied<-dummied[rownames(dist.m),-1]
    co.inc<-dummied%*%t(dummied) #Increment for co-occurence
    co<-co+co.inc
    
    #add<-data.frame(iter=rep(i, size), ind=names(ct), cluster=ct)
    #clus.bs<-rbind(clus.bs, add)
  }
  
  
  #Create the distance matrix
  dm<-co/cs
  disim<-as.dist(1-dm)
  
  #Work out which observations are promiscuous
  dm.rep<-dm
  diag(dm.rep)<-0
  most<-apply(dm.rep,1,var)
  include<-most>quantile(most, pct.exc)
  
  #Clustering on the full data set
  
  h<-hclust(disim, method)
  
  #Remove unstable obs from the distance matrix
  
  dm<-co[include,include]/cs[include,include]
  row.names(dm)<-rownames(co)[include]
  disim<-as.dist(1-dm)
  
  #Run hierarchical clustering
  h.co<-hclust(disim, method)
  
  eval<-data.frame(obs.ind=1:nd, obs.names=names(most), max.co = most, exclude=most<=quantile(most, pct.exc))
  clus<-data.frame(obs.names=h$labels, cluster = cutree(h, k))
  m<-merge(clus, eval, by="obs.names", all.x=TRUE)
  m$cluster[m$exclude==TRUE]<-0
  
  buster<-list(bhclust=h, bhclust.ex=h.co, hclust=hclust(dist), obs.eval=m)
  
  class(buster)<-"buster"
  
  return(buster)
  
}

