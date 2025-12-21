# PhyInformR
Code and Supplemental Electronic Data files for the in review BMC Evol Biol Applications Note:

PhyInformR: analysis of phylogenetic informativeness in R

ALEX DORNBURG1*, J. Nick Fisk2, JULES TAMAGNAN3, AND JEFFREY P. TOWNSEND2,4,5

1  North Carolina Museum of Natural Sciences, Raleigh, North Carolina 27601
2 Department of Biostatistics, Yale University, New Haven, Connecticut 06510, USA
3 Center for Infectious Disease Modeling and Analysis, Yale School of Public Health, Yale University, New Haven, Connecticut 06510, USA
4 Department of Ecology and Evolutionary Biology, Yale University, New Haven, Connecticut 06525, USA
5 Program in Computational Biology and Bioinformatics, Yale University, New Haven, Connecticut 06511, USA

Electronic Data files prepared by:
Alex Dornburg
alex.dornburg@naturalsciences.org

###

Full documentation for PhyInformR is available in the manual and on the web:

carolinafishes.github.io/software/phyinformR

---

## This repository: updated / extended version

This repo keeps the original PhyInformR functions and adds a set of **modernized, interactive HTML outputs** and a few **exact performance optimizations** that do **not** change the underlying math.

If you used the old script before, you can think of this as the same core algorithms, plus:

### What’s new compared to the original (“old”) script

#### 1) Exact speedups (algorithm unchanged)
- **Rate folding (unique values + multiplicities)** is used in the new HTML workflows to avoid repeating the same computations for identical site rates.
  - This is **not an approximation**: it is mathematically identical to summing over all sites.
- A folded version of the core probability calculator is included:
  - `.fold_rates()`
  - `Approximator_folded()` (same algebra as legacy `Approximator()`)

#### 2) Better interactive plots (Plotly)
- **Trackpad/mouse-wheel zoom disabled** by default to avoid accidental mis-zoom.
- **Sticky control panels** so users can keep controls visible while scrolling.
- Cleaner hover text and “publication-first” behavior (pan instead of zoom for drag defaults).

#### 3) New interactive HTML views (experimental, actively being refined)
- **Multi-locus PI profile (HTML)**  
  `informativeness.profile_multi_html()`  
  One trace per locus, with:
  - **Top-N** visibility controls
  - **Sort modes** (peak height / peak time / closest to a chosen time)
  - **Time-of-interest** control with a draggable dashed reference line (when enabled)

- **Tree + click edge → per-locus Signal−Noise lollipop (HTML)**  
  `tree_signal_noise_multi_html()`  
  Click an internal edge midpoint on the phylogeny to update a **lollipop plot** across loci:
  - baseline at 0
  - stems up/down by **S−N = P_correct − P_wrong**
  - dot at the S−N value
  - built-in **Export SVG** button (Route A; see below)

> Note on maturity: the HTML layer is the newest part of this work. It is usable now, but it is still undergoing testing and UX polish across browsers/OSes.

#### 4) Safety checks (optional self-test)
A small optional self-test block is included to confirm that modernization overrides match baseline outputs within tolerance (when enabled).

---

## Requirements

### R packages
Core PhyInformR functions rely on standard base R plus common phylogenetics packages.

For the new interactive HTML features, you will typically need:
- `ape`
- `plotly`
- `htmlwidgets`
- `htmltools`
- `jsonlite`

Install once:
```r
install.packages(c("ape","plotly","htmlwidgets","htmltools","jsonlite"))
```

### Files
- `PhyInformR.R` — main script (includes both original and new functions in this repo)
- `PhyInformR_old.r` — reference copy of the original script

---

## Quick start

## Recommended demo data (included in this repo)

For a quick end-to-end test (rates + ultrametric tree), we recommend using the **Prumetal** example dataset that ships with this repository:

- Rates: `Data/prumetalrates.rda`
- Tree:  `Data/Prumetal_timetree.phy`

These two files are the easiest way to reproduce the interactive HTML demos (PI profiles and tree-based signal–noise views) without preparing your own inputs first.

### 1) Load functions
```r
source("PhyInformR.R")
library(ape)
```

### 2) Create a toy tree + toy rate vectors (demo data)
```r
set.seed(1)
tree <- ape::rtree(12)

# Example: 20 "loci", each with 200 site rates
rates_list <- setNames(
  lapply(1:20, function(i) rexp(200, rate = 2) ),
  sprintf("Gene_%03d", 1:20)
)
```

---

## Feature 1: Multi-locus PI profile (interactive HTML)

### Generate the widget (in RStudio Viewer or browser)
```r
w <- informativeness.profile_multi_html(
  rates_list = rates_list,
  tree       = tree,
  file       = "pi_multi.html",    # set NULL to return widget only
  selfcontained = TRUE,
  default_top_n = 10
)
w
```

### Controls you will see in the HTML
- **Sort loci by**
  - `Peak height (desc)`
  - `Peak time (latest first)`
  - `Peak closest to selected time`
- **Show top N** (slider): controls which traces are visible
- **Time of interest** (slider): moves the dashed reference line and updates the ranking (when using “closest to selected time”)

### Notes / toggles
- `default_top_n` sets the initial number of visible loci.
- `times` can be provided explicitly; otherwise branching times are used.
- The dashed reference line can be draggable if the layout enables shape dragging (implementation detail inside the function).

---

## Feature 2: Tree + click edge → Signal−Noise lollipop (interactive HTML)

This is the view that matches your “baseline at 0, stems up/down, dot at S−N” sketch.

### Generate the widget (and save to HTML)
```r
w <- tree_signal_noise_multi_html(
  rates_list = rates_list,
  tree       = tree,
  s          = 4,                  # same 's' parameter used by Approximator/psnr/pnl
  file       = "tree_sn.html",
  selfcontained = TRUE,
  default_top_n = 25
)
w
```

### How to use the HTML
- **Left panel**: phylogeny with clickable edge midpoints
- **Right panel**: lollipop plot of **S−N = P_correct − P_wrong** across loci
- Click different edges to update the lollipop plot.

### Lollipop controls
- **Sort loci by**
  - `abs(S−N) (desc)` (default)
  - `S−N (desc)`
  - `S−N (asc)`
- **Show top N**: number of loci to display in the lollipop plot
- **Export SVG**: downloads a vector graphic suitable for publication workflows

---

## Export workflow (Route A: most stable, no Python)

The HTML export feature is implemented **client-side** using Plotly, with an explicit SVG export button in the Signal−Noise view.

Recommended workflow:
1. In the HTML page, click **Export SVG**.
2. Convert to PDF using one of:
   - Open the SVG in a browser and “Print → Save as PDF” (works well on macOS/Windows)
   - Inkscape / Illustrator (if you already use them)
   - Optional R helper (below)

### Optional helper: SVG → PDF in R (optional dependency)
If you want a one-command conversion in R, you can use `rsvg`:

```r
convert_svg_to_pdf <- function(svg_file, pdf_file) {
  if (!requireNamespace("rsvg", quietly = TRUE)) {
    stop("Please install.packages('rsvg')")
  }
  rsvg::rsvg_pdf(svg_file, pdf_file)
}
```

This is optional. The core HTML export does not require Python, reticulate, or Kaleido.

---

## Notes on performance

- The new HTML signal/noise workflow computes a matrix of edge × locus values.
- Rate folding (`.fold_rates()` + `Approximator_folded()`) reduces redundant work when many sites share the same rate.
- The output is identical to the non-folded formulation (within floating tolerance), but faster on real datasets where rates repeat.

---

## Common troubleshooting

### “NodeWalker(tree) not found”
Some HTML features depend on `NodeWalker(tree)` to provide parent/daughter nodes and node times. Make sure the function is defined in your `PhyInformR.R` (this repo’s version includes it).

### HTML loads but scrolling feels broken
The updated HTML functions include page-level CSS to keep scrolling enabled and control panels sticky. If you embed the widget inside another webpage/app, that container may override scrolling.

### Accidental zoom / pan
These are intentionally minimized:
- scroll-wheel zoom is disabled in Plotly config
- modebar zoom/pan tools are removed where appropriate

---

## Citation / attribution

Please cite the original PhyInformR Applications Note (in review at the time this code was distributed) and keep the author block at the top of this README unchanged when redistributing.

---

## Contact

Original Electronic Data files prepared by:
Alex Dornburg
alex.dornburg@naturalsciences.org

Maintained and modernized updates by:
Yide Jin
jinyide0202@gmail.com OR 	yide.jin@yale.edu
