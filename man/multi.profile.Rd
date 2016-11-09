\name{multi.profile}
\alias{multi.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot phylogenetic informativeness profiles based on site rate cutoffs by rate
}
\description{
A function to plot phylogenetic informativeness profiles based on site rate breaks
}
\usage{
multi.profile(rate.vector, tree, breaks, values = "display")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rate.vector}{
site rate vector, converted to class 'matrix'
}
  \item{tree}{
A phylogenetic tree
}
  \item{breaks}{
bounds for each profile
}
  \item{values}{
'display' or 'off': whether to return PI values for each time
}
}
\references{
Townsend, J. P. “Profiling Phylogenetic Informativeness” Systematic biology 56, no. 2 (2007): 222–231.
}
\author{
Alex Dornburg
}
\examples{
library("ape")
library("splines")
read.tree(system.file("extdata","polypterus_trees.phy",package="PhyInformR"))->trees
trees[[1]]->tree
as.matrix(rag1)->rates
lower<-c(0.,0.003)
upper<-c(0.003000000001,10)
breaks<-cbind(lower,upper)
multi.profile(rates, tree, breaks, values = "off")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
