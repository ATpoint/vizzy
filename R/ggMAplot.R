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
#' @param xlab x label
#' @param ylab y label
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
#' 
#' @author Alexander Toenges
#' 
#' @details 
#' => towards xval.thresh/yval.threshold/pval.threshold
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
#' => towards quantiles.x/y
#' By default in MA-plot mode the yaxis (fold changes) is automatically scaled to avoid overly wide
#' limits due to outliers. The 0.001th and 0.999th quantile of the \code{yvals} is used to set the limits.
#' Points beyond these limits will be trimmed back to these exact limits and displayed as trianges rather
#' than dots. The user can set custom quantiles or simply use an explicit \code{xlim} or \code{ylim}
#' value to manually set axis limits. The automatic axis limit setting can be turned off by setting
#' \code{quantiles.x} or \code{quantiles.y} to NULL or \code{c(0,1)}.
#' 
#' @examples 
#' dds <- DESeq2::DESeq(DESeq2::makeExampleDESeqDataSet(5000,10))
#' res <- DESeq2::results(dds) %>% data.frame %>% na.omit
#' res$baseMean<-log2(res$baseMean+1)
#' theme_set(theme_classic(base_size = 12.5))
#' 
#' # a standard MA-plot (pvalue instead of padj simply because then we get more genes in color)
#' ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$pvalue,
#' title = "DE results", subtitle = "most basic type of MA-plot")
#' 
#' # with an abline at abs(log2(1.5)) and at 0:
#' ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$pvalue,
#' y.ablines = c(log2(1.5),-log2(1.5)), title = "DE results", subtitle = "with ablines")
#' 
#' #/ show only from -2 to 2 on the y-axis, points beyond the limits will be printed as triangles:
#' ggMAplot(xval = res$baseMean, yval = res$log2FoldChange, pval = res$pvalue,
#' y.ablines = c(log2(1.5),-log2(1.5)), title = "DE results", subtitle = "with ablines",
#' ylim=c(-2,2))
#' 
#' # as Volcano:
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
#' title = "DE results", subtitle = "Volcano plot", preset = "volcano")
#' 
#' # Volcano with only abs(logFC) > 2 colored
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
#' title = "DE results", subtitle = "Volcano plot", preset = "volcano", xval.thresh=2,
#' x.ablines=c(-2,2))
#' 
#' @export
ggMAplot  <- function(xval, yval, pval=NULL,
                      xval.thresh=NULL, yval.thresh=0, pval.thresh=0.05,
                      col.base="grey50", col.up="firebrick", col.down="darkblue",
                      Up="Up", Down="Down", NonSig="NonSig",
                      point.size=1, point.alpha=0.75,
                      xlab="baseMean", ylab="logFC",
                      xlim=NULL, ylim=NULL,
                      title=waiver(), subtitle=waiver(),
                      x.ablines=NULL, x.ablines.col="black", x.ablines.lty="dashed",
                      y.ablines=NULL, y.ablines.col="black", y.ablines.lty="dashed",
                      quantiles.x=c(0, 1), quantiles.y=c(0.001, 0.999),
                      preset=c("maplot", "volcano"), no.legend=FALSE, legend.position="bottom")
{
  
  #----------------------------------
  # checks
  
  invisible(match.arg(class(pval), c("numeric", "NULL")))
  invisible(match.arg(class(xval.thresh), c("numeric", "NULL")))
  invisible(match.arg(class(yval.thresh), c("numeric", "NULL")))
  invisible(match.arg(class(pval.thresh), c("numeric", "NULL")))
  invisible(match.arg(preset, c("maplot", "volcano")))
  invisible(match.arg(legend.position, c("left", "right", "top", "bottom")))
  preset <- match.arg(preset)
  
  #/ maplot preset is just the default but volcano is not:
  if(preset=="volcano"){
    if(is.null(xval.thresh)) xval.thresh <- 0
    if(yval.thresh==0) yval.thresh <- NULL
    if(xlab=="baseMean") xlab <- "logFC"
    if(ylab=="logFC") ylab <- "-log10(FDR)"
    if(all(quantiles.x==c(0,1))) quantiles.x <- c(0.001, 0.999)
    if(all(quantiles.y==c(0.001,0.999))) quantiles.y <- c(0,1)
  }
  
  #----------------------------------
  # axis limits depending on the
  # xval and yval thresholds and
  # quantiles
  
  if(is.null(xlim)){
    xlim <- as.numeric(quantile(xval, quantiles.x))
  }
  
  if(is.null(ylim)){
    ylim <- as.numeric(quantile(yval, quantiles.y))
  }
  
  #----------------------------------
  # mutate the data into final form
  df <-
    data.frame(xval, yval) %>%
    
    #/ add p-values
    mutate(pval=case_when(is.null(pval) ~ rep(1, length(xval)), TRUE ~ pval)) %>%
    
    #/ if beyond xlim/ylim trim to specified x/ylim and plot as triangle:
    mutate(outlier=case_when(xval > xlim[2] ~ "yes_x_top",
                             xval < xlim[1] ~ "yes_x_bottom",
                             yval > ylim[2] ~ "yes_y_top",
                             yval < ylim[1] ~ "yes_y_bottom",
                             TRUE           ~ "no")) %>%
    
    mutate(xval=case_when(outlier == "yes_x_top"    ~ xlim[2],
                          outlier == "yes_x_bottom" ~ xlim[1],
                          TRUE                      ~ xval)) %>%
    
    mutate(yval=case_when(outlier == "yes_y_top"    ~ ylim[2],
                          outlier == "yes_y_bottom" ~ ylim[1],
                          TRUE                      ~ yval))
  
  #/ classify significant points depending on preset:
  if(preset=="maplot"){
    
    if(is.null(xval.thresh)){
      xval.thresh <- 0
    }
    if(is.null(yval.thresh)){
      yval.thresh <- 0
    }
    
    df <- df %>%
      mutate(signif=factor(case_when(
      pval < pval.thresh & xval > xval.thresh & yval > yval.thresh ~ get(Up),
      pval < pval.thresh & xval > xval.thresh & yval < yval.thresh ~ get(Down),
      TRUE ~ as.character(get(NonSig))),
      levels=c(Up, Down, NonSig)))
    
  } else {
    
    if(is.null(xval.thresh)){
      xval.thresh <- 0
    }
    
    df <- df %>%
      mutate(signif=factor(case_when(
        pval < pval.thresh & xval > +xval.thresh ~ get(Up),
        pval < pval.thresh & xval < -xval.thresh ~ get(Down),
        TRUE ~ as.character(get(NonSig))),
        levels=c(Up, Down, NonSig)))
    
  }
    
  #/ make sure non-outliers are dots and all others are trianges
  df$outlier <- factor(ifelse(df$outlier=="no", "no", "yes"),
                       levels = c("no", "yes"))
  shapes <- c(20, 17)
  
  
  #----------------------------------
  # assemble plot object
  gg <- 
  ggplot(df, aes(x=xval, y=yval, color=signif, shape=outlier)) + 
    geom_point() +
    scale_color_manual(values=c(col.up, col.down, col.base), drop=FALSE) +
    scale_shape_manual(values=shapes, drop=FALSE) +
    guides(shape=FALSE) +
    ggtitle(title, subtitle) +
    xlab(xlab) + ylab(ylab) +
    xlim(xlim) + ylim(ylim)
  
  #----------------------------------
  # optional ablines
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
  
  #/ optinal remove legend
  gg <- gg + theme(legend.position=legend.position,
                   legend.title=element_blank())
  
  if(no.legend) gg <- gg + theme(legend.position="none")
  
  return(gg)
  
}
