# vizzy

vizzy is a simple little R package that provides functions for common visualization tasks
related to bioinformatics, under development. Documentation will follow.
All plots are ggplot-based.

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

A brief overview over the functions. More functions will be added soon.  

In general there is no hardcoded ggplot theme so users should set `theme_set()` before
calling plots to set the desired theme.

### ggMAplot
This is a function for MAplots and Volcanos, optionally with highlighting of significant genes based on a vector of p-values.
Minimal input is a vector with average expression values (xval) and fold changes (yval). For the color highlighting a vector with pvalues (pval) is necessary. Points with pval below a user-defined threshold (default .05) are then highlighted and the number of genes being up- or downregulated is printed in the legend. 
By default the function tries to set reasonable ylimits using quantiles of the data to avoid outliers causing excessive limits. The same can optionally be done for the xaxis. All options can be turned off. One can also plot Volcanos (xaxis=fold change, yaxis=-log10(pvalue)) when setting `preset=volcano` and providing the respective data to the xval and yval arguments.

![MAs_Volcanos](https://i.ibb.co/fkX2tzv/ggMA.png)

### plot_profiles
This is a function to produce profiles plots starting from scores of bigwig files and genomic intervals from GRanges objects. 
It is possible to either use one bigwig file and plot the scores over one or many sets of genomic intervals, or use many bigwig files and plot the scores over
one set of genomic intervals. Optionally the function calculates confidence intervals for each profile curve using bootstrapping.

![profiles1](https://i.ibb.co/0yVxyvM/profiles.png)


