#' A ggplot wrapper to produce MAplots.
#' 
#' A ggplot wrapper for MA-plots and Volcanos
#' 
#' @param xval a vector with values for the x-axis, usually average expression in MAs and 
#' fold change in Volcanos
#' @param yval a vector with values for the y-axis, usually fold change in MAs and -log10(pvalue) 
#' in Volcanos
#' @param pval a vector with pvalue-like values. Points below the \code{pval.thresh} will be 
#' highlighted
#' @param xval.thresh threshold for x-axis values, see details
#' @param yval.thresh threshold for y-axis values, see details
#' @param pval.thresh threshold for pvalues, see details
#' @param basic.col the color for the data points
#' @param col.up use this color for points classified as Up, see details
#' @param col.down use this color for points classified as Down, see details
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
#' @param loess.span if not NULL then use this value (e.g. 0.2) as span to fit a loess regression for
#' yval~xval and add to the plot
#' @param loess.color color of loess.span regression line
#' @param preset the type of plot, "maplot" by default, or "volcano" to make a volcano plot, 
#' see examples
#' @param no.legend logical, whether to remove the legend
#' 
#' @author Alexander Toenges
#' 
#' @details 
#' => towards xval.thresh/yval.threshold/pval.threshold
#' We assume a typical differential expression results table, so with a fold change, 
#' an average expression value and some kind of p-value.
#' We can trigger color-highlighting of significant genes iF a pvalue is provided to the function. 
#' For MAplots points below this pval.thresh will then be classified as "Up" if yval > yval.thresh 
#' which is usually the fold change in MA-plots, or "Down" if < yval.thresh, or "NonSig" 
#' if > pval.thresh. xval.thresh is turned off by default for MAplots but if a threshold is provided
#' then only points with xval > xval.thresh will be considered for the classification.
#' In \code{preset="volcano"} the yval.thresh is ignored as the yaxis is simply -log10(pval),
#' therefore one should use \code{pval.thresh} for filtering.
#' 
#' => towards quantiles.x/y
#' By default in MAplot mode the yaxis (fold changes) is automatically scaled to avoid overly wide
#' limits due to outliers. The 0.001th and 0.999th quantile of the yvalues is used to set the limits.
#' Points beyond these limits will be trimmed back to these limits and displayed as trianged rather
#' than dots. The user can set custom quantiles or simply use an explicit \code{xlim} or \code{ylim}
#' value to manually set axis limits. The automatic axis limit setting can be turned off by setting
#' \code{quantiles.x} or \code{quantiles.y} to NULL or \code{c(0,1)}.
#' 
#' @examples 
#' dds <- DESeq2::DESeq(DESeq2::makeExampleDESeqDataSet(5000,10))
#' res <- DESeq2::results(dds) %>% data.frame %>% na.omit %>% mutate(baseMean=log2(baseMean+1))
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
#' # as Volcano:
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
#' title = "DE results", subtitle = "Volcano plot", preset = "volcano")
#' 
#' # with modified y-limits (values beyond axis limits appear as trianges) from 0 to 2:
#' ggMAplot(xval = res$log2FoldChange, yval = -log10(res$pvalue), pval = res$pvalue,
#' title = "DE results", subtitle = "trimmed", preset = "volcano", ylim = c(0, 2))
#' 
#' @export
ggMAplot  <- function(xval, yval, pval=NULL,
                      xval.thresh=NULL, yval.thresh=0, pval.thresh=0.05,
                      basic.col="grey50", col.up="firebrick", col.down="darkblue",
                      Up="Up", Down="Down", NonSig="NonSig",
                      point.size=1, point.alpha=0.75,
                      xlab="baseMean", ylab="logFC",
                      xlim=NULL, ylim=NULL,
                      title=waiver(), subtitle=waiver(),
                      x.ablines=NULL, x.ablines.col="black", x.ablines.lty="dashed",
                      y.ablines=NULL, y.ablines.col="black", y.ablines.lty="dashed",
                      quantiles.x=c(0, 1), quantiles.y=c(0.001, 0.999),
                      loess.span=NULL, loess.color="black", loess.lty="solid", 
                      preset=c("maplot", "volcano"), no.legend=FALSE)
{
  
  ####################################
  # Checks
  ####################################
  
  invisible(match.arg(class(pval), c("numeric", "NULL")))
  invisible(match.arg(class(xval.thresh), c("numeric", "NULL")))
  invisible(match.arg(class(yval.thresh), c("numeric", "NULL")))
  invisible(match.arg(class(pval.thresh), c("numeric", "NULL")))
  invisible(match.arg(preset, c("maplot", "volcano")))
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
  
  ####################################
  # classify
  ####################################
  
  #/ Given the x/y/p cutoffs classify points as "Up", "Down", "NonSig"
  #/ If no pvals are provided set all to 1
  if(is.null(pval)) {
    pval <- rep(1, length(xval))
    rm.Legend <- TRUE
  } else rm.Legend <- FALSE

  #/ Make data.frame from the three main elements x/y/p:
  df <- cbind(data.frame(xval, yval), pval)
  
  #/ Based on the threshold classify points into up/down/nonsig groups:
  classify <- lapply(c("xval.thresh", "yval.thresh"), function(x){
    
    nm <- gsub("\\..*", "", x)
    
    # threshold not specified return empty df:
    if(is.null(get(x))) return(rep("", nrow(df)))
    
    # else:
    tmp <- rep(NA, nrow(df))
    
    #/ specific classifier for MAs
    if(preset=="maplot"){
      if(x=="xval.thresh") {
        tmp[df$pval < pval.thresh & df$yval > abs(yval.thresh) & df[,nm] >= get(x)] <- "Up"
        tmp[df$pval < pval.thresh & df$yval < abs(yval.thresh) & df[,nm] >= get(x)] <- "Down"
      }
      if(x=="yval.thresh") {
        tmp[df$pval < pval.thresh & df[,nm] >= abs(get(x))] <- "Up"
        tmp[df$pval < pval.thresh & df[,nm] <= abs(get(x))] <- "Down"
      }
    }
    
    #/ specific classifiers for Volcanos:
    if(preset=="volcano"){
      if(x=="xval.thresh") {
        tmp[df$pval < pval.thresh & df[,nm] >= get(x)] <- "Up"
        tmp[df$pval < pval.thresh & df[,nm] <= get(x)] <- "Down"
      }
      if(x=="yval.thresh") {
        if(!is.null(yval.thresh)) {
          message("yval.thresh is ignored for Volcanos as the yaxis is simply -log10(pval).")
          message("Instead use pval.thresh for filtering based on significance.")
        }
      }
    }
    
    tmp[is.na(tmp)] <- "NonSig"
    return(tmp)
    
  }) %>% do.call(cbind, .)
  
  #/ Combine the columns, if both are not the same then it is NonSig,
  #/ else it is Up or Down:
  classify <- classify[,!classify[1,] %in% "",drop=FALSE]
  if(ncol(classify)>1) classify[!classify[,1]==classify[,2],1] <- "NonSig"
  tmp<-classify[,1,drop=TRUE]
    
  df$Color=factor(tmp, levels=c("Up", "Down", "NonSig"))
  
  is.up <- length(grep("Up", tmp))
  is.down <- length(grep("Down", tmp))
  is.ns <- length(grep("NonSig", tmp))
  is.total <- nrow(df)
  
  #/ Optionally fit a loess as y~x to visualize ratio yval > / < 0
  if(!is.null(loess.span)) {
    lfit <- limma::loessFit(yval, xval, method="lowess", span=loess.span)
    df<-cbind(df, data.frame(X=xval, Y=lfit$fitted))
  }
  
  ####################################
  # Winsorize axis limits
  ####################################
  
  #/ By default the y-axis limits will be winsorized by the 0.001st and 0.999th quantile
  #/ to avoid excessive axis limits due to outliers:
  df$Outliers <- rep("no", nrow(df))
  for(i in c("x", "y")){
    if(is.null(get(paste0(i, "lim")))){
      quanty <- get(paste0("quantiles.", i))
      if(is.null(quanty)) quanty <- c(0,1)
      q.top    <- as.numeric(quantile(df[,paste0(i,"val")], quanty[2]))
      q.bottom <- as.numeric(quantile(df[,paste0(i,"val")], quanty[1]))
    } else {
      q.top <- get(paste0(i, "lim"))[2]
      q.bottom <- get(paste0(i, "lim"))[1]
    }
    w.top    <- which(df[,paste0(i, "val")] > q.top)
    w.bottom <- which(df[,paste0(i, "val")] < q.bottom)
    df$Outliers[w.top] <- "yes"
    df$Outliers[w.bottom] <- "yes"
    # with top outliers as "top" and bottom outliers as "bottom" change their y-values:
    df[w.top,paste0(i, "val")] <- q.top
    df[w.bottom,paste0(i, "val")] <- q.bottom
  }; suppressWarnings(rm(q.top,w.top,q.bottom,w.bottom))
  df$Outliers <- factor(as.character(df$Outliers), levels=c("yes", "no"))
  
  # This df2 ensures that always all factor levels for Color and Outliers will be present 
  # when creating the toplevel ggplot. The actual data for geom_point will later be read from df,
  # not from df2 but even if df contains no "is.up" and no "is.down" then still df2 ensures 
  # that the legend contains all three levels for Up, Down and NonSig,
  # Idea from adapted from https://stackoverflow.com/questions/22276761/
  l <- c("Up", "Down", "NonSig")
  df2 <- data.frame(xval=rep(0,3),
                    yval=rep(0,3),
                    pval=rep(0,3),
                    Color=factor(l, levels=l),
                    Outliers=factor(c("yes", "no", "no"), levels=c("no", "yes")))
  
  ####################################
  # Assemble basic plot
  ####################################
  
  # construct the ggplot with this df2 that contains only the factor levels and some dummy data
  gg <- ggplot(data=df2, aes(x=xval, y=yval, color=Color, shape=Outliers)) + 
          geom_blank() + xlab(xlab) + ylab(ylab) +
          geom_point(alpha=point.alpha, size=point.size, data=df) +
          scale_color_manual(values=c("Up"=col.up,
                                      "Down"=col.down,
                                      "NonSig"=basic.col),
                             labels=c(paste(Up, is.up),
                                      paste(Down, is.down),
                                      paste(NonSig, is.ns))) +
          guides(shape=FALSE) +
          theme(legend.position="bottom", legend.justification="left", 
                legend.margin=margin(0, 0, 0, 0), legend.spacing.x=unit(5, "pt"),
                legend.title=element_blank()) +
          ggtitle(label=title, subtitle=subtitle)
  
  ####################################
  # loess
  ####################################
  if(!is.null(loess.span)) gg <- gg + geom_line(data=df, aes(x=X, y=Y), color=loess.color, linetype=loess.lty)
  
  ####################################
  # ablines
  ####################################
  
  #/ optional ablines:
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
  
  ####################################
  # manual limits and legend
  ####################################
  
  #/ optional axis limits
  if(is.null(xlim)) xlim <- c(round(min(df$xval)-.5,1), round(max(df$xval)+.5,1))
  if(is.null(ylim)) ylim <- c(round(min(df$yval)-.5,1), round(max(df$yval)+.5,1))
  gg <- gg + xlim(xlim) + ylim(ylim)
    
  #/ optinal remove legend
  if(no.legend) gg <- gg + theme(legend.position="none")
  
  return(gg)
  
}