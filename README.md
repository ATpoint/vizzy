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

### ggMAplot
This is a function for MAplots and Volcanos, optionally with highlighting
of significant points based on a vector of p-values. Example:  

![MAs_Volcanos](https://i.ibb.co/fkX2tzv/ggMA.png)

### plot_profiles
This is a function to produce profiles plots starting from scores of bigwig files
and genomic intervals from GRanges objects. Optionally it calculates confidence intervals
for each curve.

![profiles1](https://i.ibb.co/0yVxyvM/profiles.png)


