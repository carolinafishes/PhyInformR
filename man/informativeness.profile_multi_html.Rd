\name{informativeness.profile_multi_html}
\alias{informativeness.profile_multi_html}
\title{Multi-locus interactive phylogenetic informativeness profiles (HTML)}
\description{
Render multiple loci/partitions as an interactive phylogenetic informativeness (PI) profile
page, with per-locus traces, hover summaries, and lightweight ranking controls.
}
\usage{
informativeness.profile_multi_html(rates_list, tree, times = NULL, file = NULL, selfcontained = TRUE, default_top_n = 50)
}
\arguments{
  \item{rates_list}{Named list of numeric rate vectors, one per locus/partition. If unnamed, loci are auto-named.}
  \item{tree}{A rooted \code{\link[ape:read.tree]{phylo}} object with branch lengths (typically ultrametric).}
  \item{times}{Optional numeric vector of time points at which PI is evaluated. If \code{NULL}, uses branching times (plus 0).}
  \item{file}{Optional path to save a standalone HTML file via \code{htmlwidgets::saveWidget}. If \code{NULL}, returns the widget object.}
  \item{selfcontained}{Logical; if \code{TRUE}, embeds dependencies in a single HTML file.}
  \item{default_top_n}{Integer; number of traces shown initially (others are hidden but available via legend).}
}
\details{
This function is designed for multi-locus or multi-partition workflows where comparing
informativeness curves is the primary goal. Loci can be toggled in the legend, and
the page can expose simple ranking/sorting behaviors (e.g., by peak height or peak time)
without recomputing curves.

\strong{Status:} HTML UI features are under active testing and may evolve (layout,
controls, and export behavior may be refined).
}
\value{
A Plotly \code{htmlwidget}. If \code{file} is provided, the widget is also saved to disk.
}
\examples{
\dontrun{
data(sample.tree)
data(sample.rates)

## Demo: split one long rate vector into multiple pseudo-loci
rates_list <- split(sample.rates, ceiling(seq_along(sample.rates) / 2000))
names(rates_list) <- sprintf("Gene_\%03d", seq_along(rates_list))

informativeness.profile_multi_html(
  rates_list = rates_list,
  tree = sample.tree,
  file = "PI_multi_locus.html",
  default_top_n = 25
)
}
}
