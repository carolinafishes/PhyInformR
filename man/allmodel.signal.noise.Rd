\name{allmodel.signal.noise}
\alias{allmodel.signal.noise}
\title{
A function to compute results of QIRP, QIHP, and QIPP accomodating uneven branch lengths and user defined substitution models or base frequency distributions
}
\description{
This function calculates QIRP, QIHP, and QIPP across an uneven quartet while allowing for a user to define a substitution model and base frequencies. 
}
\usage{
allmodel.signal.noise(a, b, c, d, e, f, internode, Pi_T, Pi_C, Pi_A, Pi_G, rate_vector)
}
\arguments{
  \item{a}{
substitution model parameter. a=rTC
}
  \item{b}{
substitution model parameter. b=rAT
}
  \item{c}{
substitution model parameter. c=rGT
}
  \item{d}{
substitution model parameter
   }
  \item{e}{
substitution model parameter. e=rGC
}
  \item{f}{
substitution model parameter. f=rGA
}

  \item{internode}{
Four quartet branch lengths and internode c(62.49,62.49,65,65,3.2)->internode
}
  \item{Pi_T}{
Base Frequency of T
}
  \item{Pi_C}{
Base Frequency of C
}
  \item{Pi_A}{
Base Frequency of A
}
  \item{Pi_G}{
Base Frequency of G
}
  \item{rate_vector}{
site rate vector, converted to class 'matrix'
}
}
\details{
Please note that base frequencies must sum to 1. 
}
\value{
Function returns a matrix of QIRP, QIHP. and QIPP values
}
\references{
Townsend, J. P., Su, Z., and Tekle, Y. I. “Phylogenetic Signal and Noise: Predicting the Power of a Data Set to Resolve Phylogeny” Systematic biology 61, no. 5 (2012): 835–849.
Su, Z., Zhuo, S., Zheng, W., Francesc, L.-G., and Townsend, J. P. “The Impact of Incorporating Molecular Evolutionary Model into Predictions of Phylogenetic Signal and Noise” Frontiers in Ecology and Evolution 2, (2014): doi:10.3389/fevo.2014.00011, Available at http://dx.doi.org/10.3389/fevo.2014.00011
Su, Z. and Townsend, J. P. “Utility of Characters Evolving at Diverse Rates of Evolution to Resolve Quartet Trees with Unequal Branch Lengths: Analytical Predictions of Long-Branch Effects” BMC evolutionary biology 15, (2015): 86.
}
\author{
A. Dornburg and J.N. Fisk
}




\examples{
as.matrix(rag1)->rates
a<-1
b<-1
c<-1
d<-1
e<-1
f<-1
Pi_T<-.25
Pi_C<-.25
Pi_A<-.25
Pi_G<-.25
internode<-c(12,12,15,15,.37)
allmodel.signal.noise(a, b, c, d, e, f, internode, Pi_T, Pi_C, Pi_A, Pi_G, rates)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% htest
