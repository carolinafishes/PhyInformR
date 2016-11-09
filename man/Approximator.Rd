\name{Approximator}
\alias{Approximator}
\title{Quantify Quartet Resolution Probabilities Using 2012 Formulation}
\description{
	Quantify QIRP, QIHP, or QIPP using the equations from Townsend et al. 2012.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Approximator(t, t0, rateVector, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Time from tip of tree to focal internode }
  \item{t0}{ Focal internode length }
  \item{rateVector}{ An object containing a vector of site rates transformed to class "matrix" }
  \item{s}{ A number representing the character state space that generated the site rates (e.g., s=2 for binary data)}
}
\references{ Townsend, J. P., Su, Z., and Tekle, Y. I. “Phylogenetic Signal and Noise: Predicting the Power of a Data Set to Resolve Phylogeny” Systematic biology 61, no. 5 (2012): 835–849. }
\author{ A. Dornburg
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

as.matrix(rag1)->rr
Approximator(100,0.5,rr,3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
