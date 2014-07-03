\name{buster}
\alias{buster}
\title{Buster}
\usage{
buster(dist, n = 100, k, size = 0.66, method = "ward", outlier.th = 0.7)
}
\arguments{
  \item{dist}{A distance object}

  \item{n}{The number of times the data should be
  resampled}

  \item{k}{The number of clusters}

  \item{size}{The percentage of the original data in each
  resample}

  \item{method}{The linkage method to be passed to hclust}

  \item{outlier.th}{The threshold below which observations
  should be discarded}
}
\value{
An an object of class buster which includes an hclust
object on the co-ocurrences, an hclust object on the
original distance measure and an evalution of the stability
of the observations for model objects specific to the type
of classifier.
}
\description{
Performs bagging on the hierarchical clusters
}
\examples{
us.dist<-dist(USArrests)
bhc<-buster(us.dist, n=250, k=5, size=0.66, method='ward', outlier.th=0.64)
plot(bhc)

#Simple test

#First look at how it picks out the borderline cases

x.1<-rnorm(50, 10, 3)
y.1<-rnorm(50, 10, 3)
x.2<-rnorm(50, 12, 3)
y.2<-rnorm(50, 15, 4)
x.3<-rnorm(50, 10, 3)
y.3<-rnorm(50, 13, 4)

test.data<-data.frame(group=rep(1:3, each=50), x=c(x.1, x.2, x.3), y=c(y.1, y.2, y.3))
names<-c(paste0("group 1: ", 1:50), paste0("group 2: ", 1:50), paste0("group 3: ", 1:50))
rownames(test.data)<-names
dist<-dist(test.data[,-1])

bhc<-buster(dist, n=200, k=3, size=0.66, method='ward', outlier.th=0.9)

plot(bhc)

#Append the assigned clusters
graph.data<-cbind(test.data, bhc$obs.eval[order(bhc$obs.eval$obs.ind, decreasing=FALSE),])
plot(graph.data$x, graph.data$y, xlim=c(0,30), ylim=c(0, 30), pch = graph.data$group, col=graph.data$cluster+1)

max.co<-bhc$obs.eval$max.co[order(bhc$obs.eval$obs.ind)]
alpha<-(max.co-min(max.co))/(max(max.co)-min(max.co))
cols <- hsv(0,0,0,alpha)
plot(graph.data$x, graph.data$y, xlim=c(0,30), ylim=c(0, 30), pch = 19, col=cols)
text(graph.data$x, graph.data$y, alpha)
}
\author{
Simon Raper
}

