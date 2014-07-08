#'
#' Plot bagged hierarchical clustering
#'
#' Plot dendrograms of the original clustering and the bagged clustering
#'
#' @param buster An object of the class buster
#' @author Simon Raper
#' @examples
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
#' bhc<-buster(dist, n=200, k=3, size=0.66, method='ward', pct.exc=0.1)
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

plot.buster<-function(buster){
  
  plot(buster$bhclust, main="Dendrogram excluding unstable observations")
  dend <- as.dendrogram(buster$hclust)
  alpha<-1-buster$obs.eval$exclude[order(buster$obs.eval$obs.ind)]
  alpha[alpha==0]<-0.2
  labels_colors(dend) <- hsv(0,0,0,alpha[buster$hclust$order])
  plot(dend, main="Dendrogram including unstable observations")
  
}

