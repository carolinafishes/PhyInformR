\name{informativeness.profile}
\alias{informativeness.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot phylogenetic informativeness profiles
}
\description{
Plot phylogenetic informativeness profiles
}
\usage{
informativeness.profile(rate.vector, tree, codon = "FALSE", values = "display")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rate.vector}{
site rate vector, converted to class 'matrix'
}
  \item{tree}{
A phylogenetic tree
}
  \item{codon}{
True or false: whether to partition a given site rate matrix by codon
}
  \item{values}{
'display' or 'off': whether to return PI values for each time
}
}
\value{
Graphical display of phylogenetic informativeness
}
\references{
Townsend, J. P. “Profiling Phylogenetic Informativeness” Systematic biology 56, no. 2 (2007): 222–231.
}
\author{
A. Dornburg
}



\examples{
library("ape")
library("splines")
read.tree(system.file("extdata","polypterus_trees.phy",package="PhyInformR"))->trees
trees[[1]]->tree
as.matrix(rag1)->rates
informativeness.profile(rates, tree, codon = "FALSE", values = "display")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
