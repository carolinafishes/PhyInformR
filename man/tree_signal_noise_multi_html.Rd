\name{tree_signal_noise_multi_html}
\alias{tree_signal_noise_multi_html}
\title{Interactive tree + per-edge multi-locus signal-minus-noise lollipop plot (HTML)}
\description{
Create an interactive HTML page showing a phylogeny alongside a per-locus
\eqn{S-N = P_\{correct\} - P_\{wrong\}} lollipop plot. Clicking an edge on the tree updates the lollipop plot
and ranks loci according to the selected sorting option.
}
\usage{
tree_signal_noise_multi_html(rates_list, tree, s = 4, file = NULL, selfcontained = TRUE, default_top_n = 25)
}
\arguments{
  \item{rates_list}{Named list of numeric rate vectors, one per locus/partition. If unnamed, loci are auto-named.}
  \item{tree}{A rooted \code{\link[ape:read.tree]{phylo}} object with branch lengths (typically ultrametric).}
  \item{s}{Integer; the \code{s} parameter used by the legacy signal/noise approximator.}
  \item{file}{Optional path to save a standalone HTML file via \code{htmlwidgets::saveWidget}. If \code{NULL}, returns the widget object.}
  \item{selfcontained}{Logical; if \code{TRUE}, embeds dependencies in a single HTML file.}
  \item{default_top_n}{Integer; number of loci shown in the lollipop plot by default (top-ranked).}
}
\details{
\strong{Interpretation 2A (multi-locus):}
For each internal edge (internode interval) and each locus, the function computes
\eqn{P_\{correct\}}, \eqn{P_\{wrong\}}, and their difference \eqn{S-N}. The HTML output shows:
\enumerate{
  \item A clickable phylogeny (edge midpoints).
  \item A lollipop plot (baseline at 0; stems up/down; dot at \eqn{S-N}) across loci for the selected edge.
  \item A simple export control (route A): export SVG/PNG in-browser (vector SVG recommended for publication workflows).
}
The numeric computations follow the legacy approximator mathematics; performance is improved
via exact rate folding (unique rate values plus multiplicities), which does \emph{not} change results.
}
\value{
A Plotly \code{htmlwidget}. If \code{file} is provided, the widget is also saved to disk.
}
\seealso{
\code{\link{informativeness.profile_multi_html}}
}
\examples{
\dontrun{
data(sample.tree)
data(sample.rates)

## Demo: split rates into pseudo-loci for multi-locus interaction
rates_list <- split(sample.rates, ceiling(seq_along(sample.rates) / 2000))
names(rates_list) <- sprintf("Gene_\%03d", seq_along(rates_list))

tree_signal_noise_multi_html(
  rates_list = rates_list,
  tree = sample.tree,
  s = 4,
  file = "tree_signal_noise.html",
  default_top_n = 25
)
}
}
