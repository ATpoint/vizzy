# vizzy

vizzy is a simple little R package that provides functions for common visualization tasks
related to bioinformatics, under development. All plots are ggplot-based.
Type `?<function_name>` to see the help incl. example data and commands.

## Installation

```r

#/ Bioc dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install()
}
biocs <- c("genomation", "rtracklayer")
BiocManager::install(biocs)

#/ then install the package from Github:
install.packages("remotes")
remotes::install_github("ATpoint/vizzy")

```

## Functions

A brief overview over the functions. In general there is no hardcoded ggplot theme so users should set `theme_set()` before
calling plots to set the desired theme/layout, e.g. `theme_set(theme_bw(base_size=12))`.

### ggMAplot
This is a function to produce MA-plots and Volcanos, optionally with highlighting of significant genes based on a provided vector of pvalues. Minimal input for MA-plots is a vector with average expression values (xval) and log fold-changes (yval). For Volcanos it would be log fold-changes (xval) and `-log10(pvalue)`. For the color highlighting a vector with pvalues (pval) is necessary. Points with pval below a user-defined threshold (default .05) are then highlighted and the number of genes being up- or downregulated is printed to the legend. The user can also specify options for xval and yval thresholds to customize the highlighting. By default the function tries to set reasonable axis limits using quantiles of the data range to avoid outliers causing excessive limits, e.g. all yvals being between -2 and 2 with one outlier at -10. The user can also change axis limits manually. Points that would fall outside the limits are then trimmed and shown at exactly the chosen axis limits highlighted as trianges rather than points. This was inspired by the `DESeq2::plotMA` function. For MA-plots there is a preset `preset="maplot"` and for Volcanos there is `preset=volcano` which then. 

=> See `?ggMAplot` for a more detailed explanation. 

Here some example plots:

**MA-plots:**  
<br>

<embed src="./plots/maplots_examples.pdf" type="application/pdf">

<br>

<br>

<embed src="./plots/volcano_examples.pdf.pdf" type="application/pdf">

<br>

### plot_profiles
This is a function to produce profiles plots starting from scores of bigwig files and genomic intervals from GRanges objects. 
It is possible to either use one bigwig file and plot the scores over one or many sets of genomic intervals, or use many bigwig files and plot the scores over one set of genomic intervals. Optionally the function calculates confidence intervals for each profile curve using bootstrapping.

![profiles1](https://i.ibb.co/0yVxyvM/profiles.png)


