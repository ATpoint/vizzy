#' A ggplot wrapper to produce MAplots.
#' 
#' A ggplot wrapper for MA-plots and Volcanos
#' 
#' @param xval a vector with values for the x-axis, usually average expression in MAs and 
#' log fold-change in Volcanos
#' @param yval a vector with values for the y-axis, usually log fold-change in MAs and -log10(padj) 
#' in Volcanos
#' @param pval a vector with pvalues-like values. Points below the \code{pval.thresh} will be 
#' highlighted with color
#' @param labels a vector of equal length containing per-point labels (e.g. the gene names) which then can
#' be used downstream to label the points, see examples and details
#' @param xval.thresh absolute threshold value for x-axis values to be color-highlighted, see details
#' @param yval.thresh absolute threshold value for y-axis values to be color-highlighted, see details
#' @param pval.thresh threshold for pvalues, see details
#' @param col.base the base color for the data points, will also be used to color non-significant points
#' @param col.up use this color for points classified as "Up", see details
#' @param col.down use this color for points classified as "Down", see details
#' @param Up legend label for points being "Up"
#' @param Down legend label for points being "Down"
#' @param NonSig legend label for points being "NonSig"
#' @param point.size size of data points
#' @param point.alpha point opacity
#' @param xlab x-axis label, if NULl then defaults to "baseMean" for MA-plots and "log2FoldChange" for Volcanos
#' @param ylab y-axis label, if NULl then defaults to "log2FoldChange" for MA-plots and "-log10(padj)" for Volcanos
#' @param xlim xaxis limits, by default full data range
#' @param ylim yaxis limits, by default the 0.001th and 0.999th quantile of the data range, see details
#' @param title main title text
#' @param subtitle subtitle text
#' @param x.ablines a vector with numeric values, will use to draw vertical ablines
#' @param x.ablines.col color of x.ablines
#' @param x.ablines.lty line type of x.ablines
#' @param y.ablines a vector with numeric values, will use to draw horizontal ablines
#' @param y.ablines.col color of y.ablines
#' @param y.ablines.lty line type of y.ablines
#' @param quantiles.x use these quantiles of the xaxis data range to determine automatic axis limits, 
#' see details
#' @param quantiles.y use these quantiles of the yaxis data range to determine automatic axis limits, 
#' see details
#' @param preset the type of plot, "maplot" by default, or "volcano" to make a volcano plot, 
#' see examples
#' @param no.legend logical, whether to remove the legend
#' @param legend.position legend positon, one of left/right/top/bottom 
#' @param legend.addnumbers logical, whether to add the number of Up/Down/NonSig genes to the legend
#' @param return.data logical, whether to return the long-form data.frame that was sed for plotting.
#' In this case the output is a list with elements "plot" and "data".
#' 
#' @author Alexander Toenges
#' 
#' @details 
#' => towards xval.thresh/yval.threshold/pval.threshold \cr 
#' We assume a typical differential expression results table, so with log fold-changes, 
#' average expression values (baseMean) and some kind of p-values.
#' We can trigger color-highlighting of significant genes if a p-value is provided to the function. 
#' For MA-plots points below this pval.thresh will then be classified as "Up" if yval > yval.thresh 
#' which is usually the fold change in MA-plots, or "Down" if < yval.thresh, or "NonSig" 
#' if > pval.thresh. xval.thresh is turned off by default for MA-plots but if a threshold is provided
#' then only points with xval > xval.thresh will be considered for the classification.
#' In \code{preset="volcano"} the yval.thresh is ignored as the yaxis is simply -log10(pval),
#' therefore one should use \code{pval.thresh} for filtering.
#' 
#' => towards quantiles.x/y \cr 
#' By default in MA-plot mode the yaxis (fold changes) is automatically scaled to avoid overly wide
#' limits due to outliers. The 0.001th and 0.999th quantile of the \code{yvals} is used to set the limits.
#' Points beyond these limits will be trimmed back to these exact limits and displayed as trianges rather
#' than dots. The user can set custom quantiles or simply use an explicit \code{xlim} or \code{ylim}
#' value to manually set axis limits. The automatic axis limit setting can be turned off by setting
#' \code{quantiles.x} or \code{quantiles.y} t0 \code{c(0,1)}. If \code{xlim} or \code{ylim} are provided then
#' any quantile parameters will be overwritten.
#' 
#' => towards labels \cr 
#' The user can provide a vector of gene names via \code{labels} which is then internally added to the ggplot
#' object. That allows to manually label data points e.g. with \code{ggrepel}. We do not have an 
#' in-built option for labeling points in order to avoid overly many options, and (based on our experience)
#' adding labels to plots with many data points anyway requires custom tweaking to produce a "good-looking" plot.
#' It can be done via e.g. \code{ggrepel::geom_text_repel()} without much effort though, see last example.
#' Basically, one uses an ifelse statement that returns the gene name if it is to be plotted or an empty string if not.
#' The label is then added to the data.frame that is passed to ggplot, so one simply has to add "+ geom_label_repel" in order
#' to plot the label. One can then add ggrepel options at will, e.g. nudge_x/y to increase connector length.
#' 
#' @examples 
#' set.seed(1)
#' dds <- DESeq2::DESeq(DESeq2::makeExampleDESeqDataSet(5000,10))
#' res <- DESeq2::results(dds) %>% data.frame %>% na.omit
#' res$baseMean<-log2(res$baseMean+1)
#' res$padj <- res$pvalue
#' theme_set(theme_classic(base_size = 12.5))
#' 
#' # a standard MA-plot 
#' ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$padj,
#' title = "MA-plot", subtitle = "most basic version")
#' 
#' # with an abline at abs(log2(1.5)) and at 0:
#' ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$padj,
#' y.ablines = c(log2(1.5),-log2(1.5)), title = "MA-plot", subtitle = "with some ablines")
#' 
#' #/ show only from -2 to 2 on the y-axis, points beyond the limits will be printed as triangles:
#' ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$padj,
#' y.ablines = c(log2(1.5),-log2(1.5)), title = "MA-plot", subtitle = "with ablines and restricted at the y-axis to -2/+2",
#' ylim=c(-2,2))
#' 
#' # as Volcano:
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$padj), pval = res$padj,
#' title = "Volcano", subtitle = "most basic version", preset = "volcano")
#' 
#' # Volcano with only abs(logFC) > 2 colored
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$padj), pval = res$padj,
#' title = "Volcano", subtitle = "with some ablines", preset = "volcano", xval.thresh=2,
#' x.ablines=c(-2,2))
#' 
#' # Label significant regions, threshold padj<0.001 & abs(log2FoldChange) > log2(2).
#' # Note that for the \code{labels} option we make a custom vector, same order as for the
#' # xval/yval/pval vectors, setting those genes we do **not** want to plot as \code{""}.
#' # It is now on the user to tweak the number of genes to label and to play with \code{max.overlaps}
#' # in \code{geom_text_repel} until the plot looks satisfying.
#' library(ggrepel)
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$padj), pval = res$padj,
#' labels=ifelse(res$padj < 0.001 & abs(res$log2FoldChange) > log2(2), rownames(res), ""),
#' title = "Volcano", subtitle = "with padj < 0.001 labelled", preset = "volcano", xval.thresh=log2(2)) +
#' geom_label_repel(show.legend=FALSE, max.overlaps=20, min.segment.length=0, nudge_x=1, nudge_y=-1)
#' 
#' # or only gene number 1, see also https://ggrepel.slowkow.com/articles/examples.html for many options to customize ggrepel:
#' library(ggrepel)
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$padj), pval = res$padj,
#' labels=ifelse(rownames(res) == "gene1", rownames(res), ""),
#' title = "Volcano", subtitle = "with specific genes labelled", preset = "volcano", xval.thresh=log2(2)) +
#' geom_label_repel(aes(label=labels), show.legend=FALSE, min.segment.length=0, nudge_x=-1)
#' 
#' 
#' @export
ggMAplot  <- function(xval, yval, pval=NULL, labels=NULL,
                      xval.thresh=NULL, yval.thresh=0, pval.thresh=0.05,
                      col.base="grey50", col.up="firebrick", col.down="darkblue",
                      Up="Up", Down="Down", NonSig="NonSig",
                      point.size=2, point.alpha=0.75,
                      xlab=NULL, ylab=NULL,
                      xlim=NULL, ylim=NULL,
                      title=waiver(), subtitle=waiver(),
                      x.ablines=NULL, x.ablines.col="black", x.ablines.lty="dashed",
                      y.ablines=NULL, y.ablines.col="black", y.ablines.lty="dashed",
                      quantiles.x=c(0, 1),
                      quantiles.y=c(.001, .999),
                      preset=c("maplot", "volcano"),
                      no.legend=FALSE, legend.position="bottom", legend.addnumbers=TRUE,
                      return.data=FALSE)
{
  
  #----------------------------------
  # checks
  
  invisible(match.arg(class(pval), c("numeric", "NULL")))
  invisible(match.arg(class(xval.thresh), c("numeric", "NULL")))
  invisible(match.arg(class(yval.thresh), c("numeric", "NULL")))
  invisible(match.arg(class(pval.thresh), c("numeric", "NULL")))
  invisible(match.arg(legend.position, c("left", "right", "top", "bottom")))
  invisible(match.arg(preset, c("maplot", "volcano"))); preset <- match.arg(preset)

  #----------------------------------
  #/ if no labels add NAs:
  if(!is.null(labels)) {
    if(length(labels)!=length(xval)){
      stop("[Error] Length of label must be the same as of xval/yval/pval")
    }
  } else labels <- rep(NA, length(xval))
  
  #----------------------------------
  # axis limits depending on the xval/yval thresholds and quantiles
  if(is.null(xlim)){
    if(all(quantiles.x == c(0,1))) {
      xlim <- c(min(xval), max(xval))
    } else xlim <- as.numeric(quantile(xval, quantiles.x))
  }
  
  if(is.null(ylim)){
    if(all(quantiles.y == c(0,1))) {
      ylim <- c(min(yval), max(yval))
    } else ylim <- as.numeric(quantile(yval, quantiles.y))
  }
  
  #----------------------------------
  # mutate the data into final long form
  df <- data.frame(xval, yval, labels) %>%
    
    #/ add p-values
    mutate(pval=case_when(is.null(pval) ~ rep(1, length(xval)), TRUE ~ pval)) %>%
    
    #/ if beyond xlim/ylim trim to specified x/ylim and plot as triangle:
    mutate(out1=case_when(xval > xlim[2] ~ "yes_x_top",
                          xval < xlim[1] ~ "yes_x_bottom",
                          TRUE ~ "no"),
           out2=case_when(yval > ylim[2] ~ "yes_y_top",
                          yval < ylim[1] ~ "yes_y_bottom",
                          TRUE ~ "no")) %>%
    
    mutate(xval=case_when(out1 == "yes_x_top"    ~ xlim[2],
                          out1 == "yes_x_bottom" ~ xlim[1],
                          TRUE                      ~ xval)) %>%
    
    mutate(yval=case_when(out2 == "yes_y_top"    ~ ylim[2],
                          out2 == "yes_y_bottom" ~ ylim[1],
                          TRUE                      ~ yval)) %>%
    
    mutate(outlier=factor(case_when(out1=="no" & out2=="no" ~ "no",
                             TRUE ~ "yes"),
                          levels = c("no", "yes"))) %>%
    
    dplyr::select(-c(out1, out2))
  
  
  #/ classify significant points depending on preset:
  if(preset=="maplot"){
    
    if(is.null(xval.thresh)){
      xval.thresh <- 0
    }
    if(is.null(yval.thresh)){
      yval.thresh <- 0
    }
    
    if(is.null(xlab)) xlab <- "baseMean"
    if(is.null(ylab)) ylab <- "log2FoldChange"
    
    df <- df %>%
      mutate(signif=factor(case_when(
      pval < pval.thresh & xval > xval.thresh & yval > abs(yval.thresh) ~ get(Up),
      pval < pval.thresh & xval > xval.thresh & yval < -abs(yval.thresh) ~ get(Down),
      TRUE ~ as.character(get(NonSig))),
      levels=c(Up, Down, NonSig)))
    
  } else {
    
    if(is.null(xval.thresh)){
      xval.thresh <- 0
    }
    
    if(is.null(xlab)) xlab <- "log2FoldChange"
    if(is.null(ylab)) ylab <- "-log10(padj)"
    
    df <- df %>%
      mutate(signif=factor(case_when(
        pval < pval.thresh & xval > abs(xval.thresh) ~ get(Up),
        pval < pval.thresh & xval < -abs(xval.thresh) ~ get(Down),
        TRUE ~ as.character(get(NonSig))),
        levels=c(Up, Down, NonSig)))
    
  }
    
  #----------------------------------
  # assemble plot object
  #/ optionally add number of DEGs to legend:
  if(legend.addnumbers){
    
    legend.labs <- 
      c(paste(Up, sum(df$signif==Up)),
        paste(Down, sum(df$signif==Down)),
        paste(NonSig, sum(df$signif==NonSig)))

  } else legend.labs <- c(Up, Down, NonSig)
  
  gg <- 
  ggplot(df, aes(x=xval, y=yval, color=signif, shape=outlier, label=labels)) + 
    geom_point(size=point.size, alpha=point.alpha) +
    scale_color_manual(values=c(col.up, col.down, col.base), 
                       labels=legend.labs,
                       drop=FALSE) +
    scale_shape_manual(values=c(20, 17), drop=FALSE) +
    guides(shape=FALSE) +
    ggtitle(title, subtitle) +
    xlab(xlab) + ylab(ylab) +
    xlim(xlim) + ylim(ylim) +
    theme(legend.position=legend.position, 
          legend.justification="left",
          legend.title=element_blank())
  
  #----------------------------------
  # optional ablines for x/y:
  for(i in c("x","y")){
    
    abl <- get(paste0(i,".ablines"))
    if(is.null(abl)) next
    
    for(a in abl){
      
      if(i=="x") {
        gg<-gg+geom_vline(xintercept=abl,
                          colour=get(paste0(i, ".ablines.col")), 
                          linetype=get(paste0(i, ".ablines.lty")))
      }
      if(i=="y") {
        gg<-gg+geom_hline(yintercept=abl,
                          colour=get(paste0(i, ".ablines.col")), 
                          linetype=get(paste0(i, ".ablines.lty")))
      }
    }
  }

  #----------------------------------
  #/ if no legend:
  if(no.legend) gg <- gg + theme(legend.position="none")
  
  if(!return.data) {
    return(gg) 
  } else return(list(plot=gg, data=df))
  
}
