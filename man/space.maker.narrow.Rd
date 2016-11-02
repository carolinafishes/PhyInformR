\name{space.maker.narrow}
\alias{space.maker.narrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Core Function required to assemble more fine-scale heatmaps of QIRP as in Prum et al. 2015}
\description{
Core function required to create heatmaps as in Prum et al. 2015.
}
\usage{
space.maker.narrow(rateVector, t, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rateVector}{
An object containing a vector of site rates transformed to class "matrix" 

}
  \item{t}{
Time of quartet (e.g, bigT from Townsend et al. 2012)
}
  \item{s}{
A number representing the character state space that generated the site rates (e.g., s=2 for binary data)
}
}

\value{
Computes a range of QIRP values

}
\references{
Townsend, J. P., Su, Z., and Tekle, Y. I. “Phylogenetic Signal and Noise: Predicting the Power of a Data Set to Resolve Phylogeny” Systematic biology 61, no. 5 (2012): 835–849.
Prum, R. O., Berv, J. S., Alex, D., Field, D. J., Townsend, J. P., Lemmon, E. M., and Lemmon, A. R. “A Comprehensive Phylogeny of Birds (Aves) Using Targeted next-Generation DNA Sequencing” Nature 526, no. 7574 (2015): 569–573.
}
\author{
A. Dornburg}

\examples{
as.matrix(rag1)->rr
space.maker.narrow(rr,15,3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
