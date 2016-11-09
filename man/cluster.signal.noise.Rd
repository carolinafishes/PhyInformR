\name{cluster.signal.noise}
\alias{cluster.signal.noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantify Quartet Resolution Probabilities Using 2012 Formulation and MCMC approach
}
\description{
Quantify Quartet Resolution Probabilities (QIRP, QIHP, and QIPP) Using 2012 Formulation and MCMC approach, outputs results to file
}
\usage{
cluster.signal.noise(t, t0, rateVector, nsims, s, filename, imagename, image = "FALSE")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
Time from tip of tree to focal internode 
}
  \item{t0}{
 Focal internode length 
}
  \item{rateVector}{
An object containing a vector of site rates transformed to class "matrix" 
}
  \item{nsims}{
number of simulations
}
  \item{s}{
A number representing the character state space that generated the site rates (e.g., s=2 for binary data)
}
  \item{filename}{
Output file name for results of calculation
}
  \item{imagename}{
Output file name for results of visualization
}
  \item{image}{
True or False if image should visualized
}
}
\details{
Please note that the two files must have different names
}
\value{
returns QIRP, QIHP, and QIPP values
}
\references{
Townsend, J. P., Su, Z., and Tekle, Y. I. “Phylogenetic Signal and Noise: Predicting the Power of a Data Set to Resolve Phylogeny” Systematic biology 61, no. 5 (2012): 835–849. }
\author{
A. Dornburg
}

\examples{
as.matrix(rag1)->rr
Approximator(100,0.5,rr,3)

cluster.signal.noise(100, 2, rr, 1000,
3, filename="output_values", 
imagename="output_graphic", image = "TRUE") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
