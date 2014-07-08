\name{buster}
\alias{buster}
\title{Buster}
\usage{
buster(dist, n = 100, k, size = 0.66, method = "ward", pct.exc = 0.1,
  low.mem = FALSE)
}
\arguments{
  \item{dist}{A distance object}

  \item{n}{The number of times the data should be
  resampled}

  \item{k}{The number of clusters}

  \item{size}{The percentage of the original data in each
  resample}

  \item{method}{The linkage method to be passed to hclust}

  \item{pct.exc}{Setting this to x excludes or highlights
  the top x % observations ranked by instability. For
  example if we set this to 0.1 then the top decile of
  unstable observations are excluded or highlighted.}

  \item{low.mem}{Setting this to true uses a slower but
  less memory intensive way of calculating the
  co-occurrences}
}
\value{
An an object of class buster which includes an hclust
object on the co-ocurrences, an hclust object on the
original distance measure and an evalution of the stability
of the observations for model objects specific to the type
of classifier.
}
\description{
Performs bagging on hierarchical clusters
}
\examples{
#Testing on the iris data set
iris.dist<-dist(iris[,1:4])
bhc<-buster(iris.dist, n=250, k=3, size=0.66, method='ward', pct.exc=0.1)
plot(bhc)

#We see the unstable observations in pink.

cluster<-bhc$obs.eval$cluster[order(bhc$obs.eval$obs.ind)]
plot(iris[,1:4], col=6-cluster, pch = rep(15:17, each=50))

#Another simple test

x.1<-rnorm(50, 10, 3)
y.1<-rnorm(50, 10, 3)
x.2<-rnorm(50, 20, 3)
y.2<-rnorm(50, 10, 3)
x.3<-rnorm(50, 13, 3)
y.3<-rnorm(50, 20, 3)

test.data<-data.frame(group=rep(1:3, each=50), x=c(x.1, x.2, x.3), y=c(y.1, y.2, y.3))
names<-c(paste0("group 1: ", 1:50), paste0("group 2: ", 1:50), paste0("group 3: ", 1:50))
rownames(test.data)<-names
dist<-dist(test.data[,-1])

bhc<-buster(dist, n=200, k=3, size=0.66, method='ward', pct.exc=0.1)

plot(bhc)

#Append the assigned clusters
graph.data<-cbind(test.data, bhc$obs.eval[order(bhc$obs.eval$obs.ind, decreasing=FALSE),])
plot(graph.data$x, graph.data$y, xlim=c(0,30), ylim=c(0, 30), pch = graph.data$group, col=graph.data$cluster+1)

max.co<-bhc$obs.eval$max.co[order(bhc$obs.eval$obs.ind)]
alpha<-(max.co-min(max.co))/(max(max.co)-min(max.co))
cols <- hsv(0,0,0,alpha)
plot(graph.data$x, graph.data$y, xlim=c(0,30), ylim=c(0, 30), pch = 19, col=cols)
}
\author{
Simon Raper
}

