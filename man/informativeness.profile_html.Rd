\name{informativeness.profile_html}
\alias{informativeness.profile_html}
\title{Interactive phylogenetic informativeness profile (HTML)}
\description{
Compute a phylogenetic informativeness (PI) profile and render it as an
interactive HTML graphic (Plotly).
}
\usage{
informativeness.profile_html(rate.vector, tree, codon = "FALSE", file = NULL, selfcontained = TRUE)
}
\arguments{
  \item{rate.vector}{Numeric vector of site rates for a locus/partition.}
  \item{tree}{A rooted \code{\link[ape:read.tree]{phylo}} object with branch lengths (typically ultrametric).}
  \item{codon}{Character flag kept for compatibility with the legacy interface
  (e.g., \code{"FALSE"}).}
  \item{file}{Optional path to save a standalone HTML file via \code{htmlwidgets::saveWidget}. If \code{NULL}, returns the widget object.}
  \item{selfcontained}{Logical; if \code{TRUE}, embeds dependencies in a single HTML file.}
}
\details{
This function is intended as a modern HTML alternative to the legacy static plotting
workflow. It uses the same PI calculations as the underlying profile functions, but
renders output with Plotly to support hover tooltips and interactive exploration.

\strong{Note on saving:} \code{selfcontained = TRUE} may require a working Pandoc
installation depending on the \code{htmlwidgets} setup on your system.
}
\value{
A Plotly \code{htmlwidget}. If \code{file} is provided, the widget is also saved to disk.
}
\examples{
\dontrun{
## Example (assuming the package provides sample data objects):
data(sample.rates)
data(sample.tree)

informativeness.profile_html(
  rate.vector = sample.rates,
  tree = sample.tree,
  codon = "FALSE",
  file = "PI_profile.html"
)
}
}
