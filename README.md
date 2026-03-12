# PhyInformR

**PhyInformR** is an R toolkit for computing and visualizing **phylogenetic informativeness (PI)**.
It helps researchers choose which genomic loci (genes, partitions) are most useful for
resolving specific parts of a phylogenetic tree, and at what time depths.

This branch (`html-testing`) adds **interactive HTML visualizations** powered by Plotly,
significant **performance improvements**, and a built-in **self-test** — while keeping
every original function fully intact.

---

## Table of Contents

1. [Background](#background-what-is-phylogenetic-informativeness)
2. [What's New in This Version](#whats-new-in-this-version)
3. [Installation](#installation)
4. [Data Included](#data-included)
5. [Quick Start](#quick-start)
6. [Function Reference](#function-reference)
   - [PI Profiles — static plots](#pi-profiles--static-plots)
   - [PI Profiles — interactive HTML](#pi-profiles--interactive-html)
   - [Multi-Locus Interactive Profile](#multi-locus-interactive-profile)
   - [Signal / Noise Analysis](#signal--noise-analysis)
   - [Bayesian Quartet Probabilities](#bayesian-quartet-probabilities)
   - [Tree Signal-Informativeness Plots](#tree-signal-informativeness-plots)
7. [Self-Test / Verification](#self-test--verification)
8. [Citation](#citation)
9. [Maintainers](#maintainers)

---

## Background: What Is Phylogenetic Informativeness?

When you build a phylogenetic tree from DNA sequences, not every site in your alignment
is equally useful.

- Sites that **never mutate** look the same in all species — they carry no information
  about relationships.
- Sites that **mutate too fast** change so many times that the original signal is
  overwritten — they can actively mislead the analysis.

**Phylogenetic informativeness (PI)** quantifies how useful a set of sites is for
resolving a *specific divergence event* at a *specific point in time* (Townsend 2007).
By plotting PI against time, you can see which genes are most informative for your
particular question — for example, resolving deep splits among bird orders vs. shallow
splits within a genus.

PhyInformR also computes **signal/noise/bias** probabilities: given a locus and a tree
topology, what is the probability of recovering the correct tree, a wrong tree, or an
unresolved polytomy?

---

## What's New in This Version

| Feature | Details |
|---------|---------|
| **Interactive HTML plots** | All major plots have an `_html` counterpart built with `plotly`. Hover for values, zoom, pan, export SVG. |
| **Multi-locus PI viewer** | `informativeness.profile_multi_html()` shows hundreds of loci at once with sort controls, top-N filter, and a draggable time-reference line. |
| **Profile-peak locus filter** | Removes loci whose PI peaks *before* the cumulative profile peak — keeping only the most relevant loci for your clade of interest. |
| **Performance overrides** | Core loops are vectorized; identical site rates are deduplicated to avoid repeating the same calculation (~10–100x faster on large datasets). |
| **Self-test mode** | Set `RUN_SELF_TEST <- TRUE` to verify that optimized code matches the original to within floating-point tolerance. |
| **No Python dependency** | `allmodel.signal.noise()` no longer calls an external Python script. Everything runs in pure R. |
| **tidyverse helpers** | `post.su()` and `su.bayes()` use `purrr` + `dplyr` instead of slow `rbind`-in-loops. |

---

## Installation

### Step 1 — Install R (if you have not already)

Download the latest version from: https://cran.r-project.org/

We also recommend **RStudio** as a beginner-friendly editor: https://posit.co/download/rstudio-desktop/

### Step 2 — Install required packages

Open R or RStudio and paste this into the Console, then press Enter:

```r
install.packages(c(
  "ape", "splines", "gplots", "RColorBrewer", "phytools",
  "foreach", "iterators", "geiger", "doParallel",
  "gridExtra", "ggplot2", "hexbin", "PBSmodelling",
  "dplyr", "purrr", "tibble",
  "plotly", "htmlwidgets", "htmltools", "jsonlite"
))
```

This may take a few minutes. Answer "yes" if asked to compile from source.

> **Tip for Windows users:** If `PBSmodelling` fails to install, try:
> `install.packages("PBSmodelling", type = "binary")`

### Step 3 — Download this script

Save `PhyInformR.R` from this repository to a folder on your computer.

### Step 4 — Load the functions

In R, set your working directory to the folder containing the script, then run:

```r
setwd("path/to/your/folder")   # change this to your actual folder path
source("PhyInformR.R")
```

All functions are now available in your R session.

---

## Data Included

The `data/` folder contains example datasets:

| File | Description |
|------|-------------|
| `prumetalrates.rda` | Site-specific substitution rates for ~250 avian loci (from Prum et al. 2015) |
| `Prumetal_timetree.phy` | Time-calibrated bird phylogeny in Newick format |

Load them like this:

```r
load("data/prumetalrates.rda")        # makes an object named prumetalrates
tree <- ape::read.tree("data/Prumetal_timetree.phy")
```

`prumetalrates` is a matrix where each **row** is one locus and each **column** is
the substitution rate for one site.

---

## Quick Start

```r
source("PhyInformR.R")

# Load the example data
load("data/prumetalrates.rda")
tree <- ape::read.tree("data/Prumetal_timetree.phy")

# 1. Static PI profile for one locus (opens a plot window)
informativeness.profile(prumetalrates[1, ], tree)

# 2. Same plot but interactive HTML (opens in your browser or saves to file)
informativeness.profile_html(prumetalrates[1, ], tree, file = "PI_locus1.html")

# 3. View many loci at once in an interactive HTML page
rates_list <- lapply(seq_len(nrow(prumetalrates)),
                     function(i) as.numeric(prumetalrates[i, ]))
names(rates_list) <- paste0("Locus_", seq_len(nrow(prumetalrates)))

informativeness.profile_multi_html(rates_list, tree, file = "PI_all_loci.html")
```

---

## Function Reference

### PI Profiles — static plots

---

#### `informativeness.profile(rate.vector, tree, codon = "FALSE", values = "display")`

Draws a two-panel figure: the phylogenetic tree on top and the PI profile curve below.

**Arguments:**

| Argument | What to provide |
|----------|----------------|
| `rate.vector` | Numeric vector of site-specific rates for one locus (one row of `prumetalrates`) |
| `tree` | A `phylo` object (from `ape::read.tree`) |
| `codon` | `"FALSE"` (default) for all sites together; `"TRUE"` to split into codon positions 1, 2, 3 |
| `values` | `"display"` returns the profile matrix; `"off"` suppresses the return value |

**Example:**
```r
informativeness.profile(prumetalrates[1, ], tree)
```

---

#### `multi.profile(rate.vector, tree, breaks, values = "display")`

Plots PI profiles for multiple **rate-based partitions** on one graph.
You define partitions by specifying rate value ranges.

**Arguments:**

| Argument | What to provide |
|----------|----------------|
| `breaks` | A matrix with 2 columns: lower bound and upper bound of each rate bin. One partition per row. |

**Example:**
```r
# Define two partitions: slow sites (rate 0–0.005) and fast sites (rate 0.005–0.02)
breaks <- matrix(c(0,    0.005,
                   0.005, 0.02), ncol = 2, byrow = TRUE)
multi.profile(prumetalrates[1, ], tree, breaks)
```

---

#### `defined.multi.profile(rate.vector, tree, breaks, values = "display")`

Same as `multi.profile`, but partitions are defined by **site index ranges**
instead of rate value ranges.

**Example:**
```r
# Partition sites 1–200 and sites 201–400
breaks <- matrix(c(1,   200,
                   201, 400), ncol = 2, byrow = TRUE)
defined.multi.profile(prumetalrates[1, ], tree, breaks)
```

---

### PI Profiles — interactive HTML

---

#### `informativeness.profile_html(rate.vector, tree, codon = "FALSE", file = NULL, selfcontained = TRUE)`

Interactive Plotly version of `informativeness.profile`.
- Hover over any point to see the exact time and PI value.
- Zoom and pan the plot.
- If `file = "output.html"`, saves a standalone HTML file you can share or open in any browser.

**Example:**
```r
informativeness.profile_html(prumetalrates[1, ], tree, file = "PI_locus1.html")
```

---

### Multi-Locus Interactive Profile

---

#### `informativeness.profile_multi_html(rates_list, tree, times = NULL, file = NULL, selfcontained = TRUE, default_top_n = 50, filter_by_profile_peak = FALSE, return_filter_info = TRUE)`

Shows PI curves for **many loci at once** in a single interactive HTML page.

**Arguments:**

| Argument | What to provide |
|----------|----------------|
| `rates_list` | A **named list** of rate vectors, one element per locus |
| `tree` | Phylogenetic tree |
| `times` | Optional custom time grid (default: all branching times from the tree) |
| `file` | Path to save the HTML file, e.g. `"PI_multi.html"` |
| `default_top_n` | Number of loci shown on first load (default: top 50 by PI peak height) |
| `filter_by_profile_peak` | If `TRUE`, hides loci whose PI peaks earlier than the summed profile peak |
| `return_filter_info` | If `TRUE`, the function returns a list reporting which loci were filtered |

**What you can do in the HTML viewer:**

- **Sort loci** by peak height, peak time, or closeness to the reference time
- **Slide the Top N control** to show more or fewer loci
- **Drag the dashed vertical line** to set the time of interest — the ranking table updates automatically
- **Click "Download SVG"** to export the figure (no Python needed)
- **Click legend entries** to show/hide individual loci

**Example:**
```r
rates_list <- lapply(seq_len(nrow(prumetalrates)),
                     function(i) as.numeric(prumetalrates[i, ]))
names(rates_list) <- paste0("Locus_", seq_len(nrow(prumetalrates)))

result <- informativeness.profile_multi_html(
  rates_list, tree,
  file = "PI_all_loci.html",
  default_top_n = 30,
  filter_by_profile_peak = TRUE   # remove loci that peak too early
)

# See which loci were filtered out:
result$removed_loci
```

---

### Signal / Noise Analysis

These functions answer the question: *"Given this locus and this tree, how likely am I to
recover the correct topology?"*

---

#### `allmodel.signal.noise(a, b, c, d, e, f, internode, Pi_T, Pi_C, Pi_A, Pi_G, rate_vector)`

Computes the probability of resolving a **four-taxon tree** correctly, incorrectly, or
as a polytomy. Uses a full **GTR substitution model**.

The four-taxon tree has this shape:

```
Leaf1 ---[ext1]---|
                  |---[internal]---|---[ext3]--- Leaf3
Leaf2 ---[ext2]---|                |---[ext4]--- Leaf4
```

**Arguments:**

| Argument | What to provide |
|----------|----------------|
| `a` | T ↔ C substitution rate |
| `b` | T ↔ A substitution rate |
| `c` | T ↔ G substitution rate |
| `d` | C ↔ A substitution rate |
| `e` | C ↔ G substitution rate |
| `f` | A ↔ G substitution rate |
| `internode` | Numeric vector of length **5**: `c(ext1, ext2, ext3, ext4, internal)` |
| `Pi_T, Pi_C, Pi_A, Pi_G` | Equilibrium base frequencies (must sum to 1.0) |
| `rate_vector` | Site-specific rates for the locus |

**Returns:** a vector `c(P_incorrect, P_polytomy, P_correct)`

**Example:**
```r
result <- allmodel.signal.noise(
  a = 1, b = 1, c = 1, d = 1, e = 1, f = 1,
  internode = c(0.1, 0.1, 0.1, 0.1, 0.05),
  Pi_T = 0.25, Pi_C = 0.25, Pi_A = 0.25, Pi_G = 0.25,
  rate_vector = as.numeric(prumetalrates[1, ])
)
cat("P(correct tree):", result[3], "\n")
```

---

#### `cluster.signal.noise(t, t0, rateVector, nsims, s, filename, imagename, image = "FALSE")`

Monte Carlo simulation of signal/noise. Runs `nsims` replicates and saves the results.

| Argument | Description |
|----------|-------------|
| `t` | Total depth of the tree (time from root to present) |
| `t0` | Time of the divergence event you are testing |
| `rateVector` | Site-specific rates |
| `nsims` | Number of simulations (larger = more accurate) |
| `s` | Number of character states (4 for DNA) |
| `filename` | Text file to save the numeric output table |
| `imagename` | PDF file to save the signal/noise bar chart |

---

#### `parallel.cluster.signal.noise(t, t0, rateVector, nsims, s, filename, imagename, image = "TRUE")`

Same as above, but uses **multiple CPU cores** for faster computation.
Set the number of cores at the top of the script:

```r
registerDoParallel(cores = 8)   # adjust to your machine
```

---

#### `graph.signal.noise_html(currentprobdist, rateVector, file = NULL, selfcontained = TRUE)`

Interactive HTML bar chart of a signal/noise probability distribution.
- **Grey bars** = probability mass for the wrong tree
- **Black bar** = polytomy
- **Blue bars** = correct tree (the signal you want)

---

### Bayesian Quartet Probabilities

These functions work with **posterior distributions of trees** from Bayesian analyses
(e.g., from MrBayes or BEAST). They calculate resolution probabilities across the
whole posterior, not just a single best tree.

---

#### `su.bayes(a, b, c, d, e, f, Pi_T, Pi_C, Pi_A, Pi_G, rate_vector, quart, tree)`

For each tree in the posterior, extracts the relevant quartet internode and computes
P(correct), P(polytomy), P(incorrect) using the GTR model.

| Argument | What to provide |
|----------|----------------|
| `quart` | Character vector of exactly 4 taxon names forming the quartet of interest |
| `tree` | A **list** of `phylo` objects (your posterior sample) |

**Returns:** a matrix — one row per posterior tree, columns for probabilities and branch lengths.

---

#### `plotPosterior(final, plotType = "QIPs")`

Visualizes the output of `su.bayes`:

- `plotType = "QIPs"` — three hexbin scatter plots of each probability vs. internode length.
  Good for seeing how resolution confidence changes with branch length.
- `plotType = "violin"` — violin + boxplot for each probability type across the posterior.

**QIRP** = Quartet Internode Resolution Probability (correct tree — you want this high)
**QIPP** = Quartet Internode Polytomy Probability (unresolved — lower is better)
**QIHP** = Quartet Internode Homoplasy Probability (wrong tree — you want this low)

---

#### `plotPosterior_html(final, plotType = "QIPs", file = NULL, selfcontained = TRUE)`

Interactive HTML version of `plotPosterior`. Hover over any point to see the exact
values for that posterior tree sample.

---

### Tree Signal-Informativeness Plots

---

#### `PlotTreeSI(tree, ratevector, s)`

Draws the phylogenetic tree with **horizontal probability segments** overlaid on each
internode. The height of each segment shows P(correct tree) for that branch given your
locus. High segments = this locus is good for resolving that internode.

---

#### `Plot.Another.TreeSI(tree, ratevector, s, color, type)`

Adds a second locus (or model) to an existing `PlotTreeSI` plot using a different
color and line type. Useful for comparing two loci side-by-side on the same tree.

**Example:**
```r
PlotTreeSI(tree, prumetalrates[1, ], s = 4)
Plot.Another.TreeSI(tree, prumetalrates[2, ], s = 4, color = "red", type = 2)
```

---

#### `space.maker(rateVector, t, s)` and `space.maker.narrow(rateVector, t, s)`

Compute a grid of correct-resolution probabilities at 20 evenly spaced divergence times
between 0 and `t` (or 0 and `t/2` for `space.maker.narrow`).
Useful for visualizing the "window of opportunity" — the time range where this locus is
informative enough to be useful.

---

## Self-Test / Verification

The script includes an internal regression test. Run it once after downloading to verify
the optimized functions match the original code to within numerical precision:

```r
# Step 1: open PhyInformR.R in a text editor and change this line near the top:
RUN_SELF_TEST <- TRUE   # was FALSE

# Step 2: source the script
source("PhyInformR.R")

# If everything passes, you will see:
# PASS: modernization overrides match baseline within tolerance.

# Step 3: set it back to FALSE for normal use
RUN_SELF_TEST <- FALSE
```

The self-test checks five core functions: `site.summer`, `inform.profile.generator2`,
`Approximator`, `ExMaker`, and `allmodel.signal.noise`.

---

## Citation

If you use PhyInformR in your research, please cite:

> Dornburg A, Fisk JN, Tamagnan J, Townsend JP (2016).
> **PhyInformR: phylogenetic experimental design and phylogenomic data exploration in R.**
> *BMC Evolutionary Biology* 16:262.
> https://doi.org/10.1186/s12862-016-0837-3

The underlying informativeness theory:

> Townsend JP (2007).
> **Profiling phylogenetic informativeness.**
> *Systematic Biology* 56(2):222–231.
> https://doi.org/10.1080/10635150701311406

The signal/noise/bias framework:

> Su Z, Townsend JP (2015).
> **Utility of characters evolving at diverse rates of evolution to resolve quartet trees
> with unequal branch lengths: analytical predictions of long-branch effects.**
> *BMC Evolutionary Biology* 15:86.
> https://doi.org/10.1186/s12862-015-0364-7

---

## Maintainers

| Role | Name | Affiliation |
|------|------|-------------|
| Original developer | Alex Dornburg | North Carolina Museum of Natural Sciences |
| Current maintainer | Yide Jin | Yale University |

For bug reports or questions, please open an Issue on GitHub.
