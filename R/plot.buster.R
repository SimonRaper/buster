#'
#' Plot bagged hierarchical clustering
#'
#' Plot dendrograms of the original clustering and the bagged clustering
#'
#' @param buster An object of the class buster
#' @author Simon Raper
#' @examples
#' 
#' #On the US arrest data
#' 
#' us.dist<-dist(USArrests)
#' bhc<-buster(us.dist, n=250, k=5, size=0.66, method='ward', outlier.th=0.9)
#' plot(bhc)
#' 
#' #Simple test
#' 
#' #First look at how it picks out the borderline cases
#' 
#' x.1<-rnorm(50, 10, 3)
#' y.1<-rnorm(50, 10, 3)
#' x.2<-rnorm(50, 20, 3)
#' y.2<-rnorm(50, 15, 4)
#' x.3<-rnorm(50, 10, 3)
#' y.3<-rnorm(50, 20, 4)
#' 
#' test.data<-data.frame(group=rep(1:3, each=50), x=c(x.1, x.2, x.3), y=c(y.1, y.2, y.3))
#' names<-c(paste0("group 1: ", 1:50), paste0("group 2: ", 1:50), paste0("group 3: ", 1:50))
#' rownames(test.data)<-names
#' dist<-dist(test.data[,-1])
#' 
#' bhc<-buster(dist, n=250, k=3, size=0.66, method='ward', outlier.th=0.68)
#' plot(bhc)

plot.buster<-function(buster){
  
  plot(buster$bhclust, main="Clustering based on co-occurrence")
  dend <- as.dendrogram(buster$hclust)
  alpha<-1-buster$obs.eval$exclude
  alpha[alpha==0]<-0.2
  labels_colors(dend) <- hsv(0,0,0,alpha[buster$hclust$order])
  plot(dend, main="Clustering based on orginal distance matrix")
  
}
