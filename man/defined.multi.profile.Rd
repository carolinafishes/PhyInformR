\name{defined.multi.profile}
\alias{defined.multi.profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to plot phylogenetic informativeness profiles based on site rate cutoffs by partition 

}
\description{
A function to plot phylogenetic informativeness profiles based on site rate cutoffs (e.g., loci)
}
\usage{
defined.multi.profile(rate.vector, tree, breaks, values = "display")
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
\value{
Graphical display of phylogenetic informativeness

}
\references{
Townsend, J. P. “Profiling Phylogenetic Informativeness” Systematic biology 56, no. 2 (2007): 222–231.
Dornburg, A., Townsend, J. P., Friedman, M., and Near, T. J. “Phylogenetic Informativeness Reconciles Ray-Finned Fish Molecular Divergence Times” BMC evolutionary biology 14, (2014): 169.
Prum, R. O., Berv, J. S., Alex, D., Field, D. J., Townsend, J. P., Lemmon, E. M., and Lemmon, A. R. “A Comprehensive Phylogeny of Birds (Aves) Using Targeted next-Generation DNA Sequencing” Nature 526, no. 7574 (2015): 569–573.
}
\author{
A. Dornburg
}
\examples{
library("ape")
library("splines")
read.tree(system.file("extdata","Prumetal_timetree.phy",package="PhyInformR"))->tree
prumetalrates->Prates
as.matrix(Prates)->rates
lower<-c(0, 1594)
upper<-c(1595,2787)
breaks<-cbind(lower,upper)
defined.multi.profile(rates, tree, breaks, values = "off")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
