\name{Plot.Another.TreeSI}
\alias{Plot.Another.TreeSI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to visualize QIRP of additional data across a tree topology
}
\description{
A function to visualize QIRP of additional data across a tree topology
}
\usage{
Plot.Another.TreeSI(tree, ratevector, s, color, type)
}
\arguments{
  \item{tree}{
A phylogenetic tree
}
  \item{ratevector}{
site rate vector, converted to class 'matrix'
}
  \item{s}{
A number representing the character state space that generated the site rates e.g., s=2 for binary data }
  \item{color}{
Color of line e.g., color="blue"
}
  \item{type}{
plotting function for lines e.g., type=3
}
}

\value{
A graphical representation of QIRP at every node

}
\references{
Hwang, J., Jonathan, H., Qi, Z., Yang, Z. L., Zheng, W., and Townsend, J. P. “Solving the Ecological Puzzle of Mycorrhizal Associations Using Data from Annotated Collections and Environmental Samples - an Example of Saddle Fungi” Environmental microbiology reports 7, no. 4 (2015): 658–667.
}
\author{
A. Dornburg
}
\examples{
library("ape")
##First plot a locus
as.matrix(rag1)->rr
read.tree(system.file("extdata","polypterus_trees.phy",package="PhyInformR"))->trees
trees[[1]]->tree
PlotTreeSI(tree,rr,3)
##Now plot another locus to compare
as.matrix(prumetalrates[1:1594])->candidate_locus
Plot.Another.TreeSI(tree,candidate_locus,3, color="red",type=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
