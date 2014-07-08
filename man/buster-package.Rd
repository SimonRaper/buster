\name{buster-package}
\alias{buster-package}
\docType{package}
\title{
Performs bagging on hierarchical clusters
}
\description{
More about what it does (maybe more than one line)
z
}
\details{
\tabular{ll}{
Package: \tab buster\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-03\cr
License: \tab What license is it under?\cr
}
y
}
\author{
Simon Raper
}
\references{
x
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
x
}
\examples{
iris.dist<-dist(iris[,1:4])
bhc<-buster(iris.dist, n=250, k=3, size=0.66, method='ward', pct.exc=0.1)
plot(bhc)
}
