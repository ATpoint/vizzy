#' Profile Plots Over Genomic Intervals
#' 
#' Produces profile plots based on scores from bigwig files and genomic intervals
#' from GRanges/GRangesList objects.
#'
#' @param bigwigs a vector with paths to bigwig files
#' @param ranges genomic intervals to use as GRanges or GRangesList, see details
#' @param names vector of same length as profiles to plot used as names for the legend, see details.
#' @param theme_set the ggplot theme, default is theme_bw
#' if one wants to use a global \code{theme_set()} then set this to NULL to turn default off.
#' @param theme_basesize the basesize if theme_set contains a valid theme
#' @param x.ticks a numeric vector that indicates the x-axis ticks to use.
#' Allowed values are 1 to \code{width(ranges)}, see details.
#' @param x.ticks.lab vector of same length as \code{x.ticks} storing the labels for each tick
#' @param x.ticks.angle numeric value, rotate tick labels by this value
#' @param x.ticks.kb logical, whether to use bp or kb as label for the default x.ticks.lab.
#' See details.
#' @param x.gridmajor logical, whether to add major vertical grid lines. Only relevant when using themes
#' that actually have grid lines, e.g. theme_bw but not theme_minimal.
#' @param center.name name of the x.ticks.lab when  x.ticks=NULL, see details.
#' @param ylab y-axis label
#' @param title main title
#' @param subtitle subtitle
#' @param colors a vector of same length as profiles to plot with valid color names.
#' If NULL then colors from a colorblind-friendly palette will be used.
#' @param loess.span numeric, the alpha parameter from loess. If not NULL then use
#' loess smoothing on each profile with this span.
#' @param legend.ncol numeric if having multiple profiles arrange legend in that many columns
#' @param lwd numeric, the line width, default is .75
#' @param CI.level numeric, if not NULL will calculate a confidence interval for every profile
#' This value is then the confidence level, default is .95 for the 95th CI.
#' @param CI.reps numeric, the number of bootstrap replicates
#' @param CI.workers workers to parallelize the bootstrapping
#' @param CI.alpha the opacity when plotting the CIs as geom_ribbon, default is .2
#' 
#' @details 
#' => There are two modes this function can use. 
#' 1) Use one bigwig file and plot over one or multiple sets of genomic intervals.
#' If using one set then ranges must be a GRanges object. If multiple then ranges must be
#' a GRangesList. If the list is named will use these names for the legend, else will simply
#' name the sets as set1, set2... etc.
#' This mode could be useful in a situation where one wants to plot e.g. one ChIP-seq sample over
#' multiple genomic sites, e.g. a set of TSS and a set of distal enhancers.
#' 2) Use many bigwig files and one set of genomic intervals. 
#' In this case \code{bigwigs} should be a vector of paths to multiple bigwigs.
#' The \code{ranges} must then be a GRanges object. If no \code{names} are provided then will
#' use the basenames of the bigwigs for the legend names.
#' 
#' => Towards the tick parameters (x.ticks etc)
#' It is important to note that the ranges in this function are always relative towards 
#' the provides genomic intervals. No matter what the genomic coordiantes are the ranges
#' will be treated as an interval from 1 to \code{width(ranges)}, so if you provide intervals
#' of width 2kb (all ranges must have the same width btw) then the default behaviour is to plot
#' these intervals using as ticks the leftpost position, the center position and the rightmost position.
#' By default these are labelled as -(width(ranges)/2), \code{center.name} and +(width(ranges)/2).
#' Means that for a 2kb interval one would get default ticks as -1kb, center, +1kb.
#' In this default mode is \code{x.ticks.kb=TRUE} (default) will label the interval in "kb" as above,
#' and if FALSE then it would be labelled as 1000bp, center, +1000bp. Option only relevant if \code{x.ticks=NULL}.
#' The \code{center.name} argument then determines the name of "center" in this above example. 
#' Not relevant if using custom ticks.
#' The ticks can be modified to any value between 1 and the width of the interval. 
#' E.g. \code{x.ticks<-c(1,500,2000)} will set the ticks to the start of the interval, 500bp and 2kb.
#' The vector can have any length. The \code{x.ticks.lab} argument can be used to provide the names for
#' these ticks, must have same length as \code{x.ticks}.
#' 
#' => Towards the confidence intervals
#' If this option is turned on then function will use bootstrapping to estimate confidence intervals
#' for every profile. Technically the profiles are the colMeans of the scores for every basepair of 
#' the genomic interval. The function will then bootstrap the scores to estimate CIs for these colMeans
#' and add it to the plot as via geom_ribbon. The bootstrapping is repeated \code{CI}.reps times. 
#' This is all done with the boot package with code adapted from the Bioconductor package ChIPseeker.
#' This does not make sense together with \code{loess.span} so with loess smoothing.
#' 
#' @author Alexander Toenges
#'  
#' @examples 
#' # an example plot
#' ranges=readRDS(paste0(system.file("extdata",package="vizzy"), "/example_ranges.rds"))
#' ranges2=GRangesList(group1=ranges[1:50],group2=ranges[51:100])
#' bigwigs=list.files(system.file("extdata", package="vizzy"), pattern = ".bigwig",
#'                    full.names = TRUE)
#'          
#' # multiple bigwigs over one set of regions with CIs     
#' plot_profiles(bigwigs=bigwigs, ranges=ranges, CI.level=.95, CI.reps=50)
#' 
#' # one bigwigs over multiple sets of peaks with CIs
#' plot_profiles(bigwigs=bigwigs[1], ranges=ranges2, CI.level=.95, CI.reps=50)
#' 
#' # custom ticks, could be positions relative to the start of the ranges being a TSS:
#' my.ticks <- c(1, 200, 500, 1000, 2000)
#' plot_profiles(bigwigs=bigwigs, ranges=ranges, 
#'               x.ticks=my.ticks, x.ticks.angle=45,
#'               x.ticks.lab=c("TSS", "+0.2kb", "+0.5kb", "+1kb", "+2kb"))
#' 
#' @references 
#' # Parsing bigwigs/GRanges into a basepair-resolution ScoreMatrix:
#' genomation: a toolkit to summarize, annotate and visualize genomic intervals
#' Akalin, A et al (2014) Bioinformatics, Volume 31, Issue 7, 1 April 2015, Pages 1127–1129
#' https://doi.org/10.1093/bioinformatics/btu775
#' 
#' # colorblind-friendly color palette:
#' dittoSeq: universal user-friendly single-cell and bulk RNA sequencing visualization toolkit
#' Bunis et al (2020) Bioinformatics, btaa1011, https://doi.org/10.1093/bioinformatics/btaa1011
#
#' @export
plot_profiles <- function(bigwigs, ranges, names=NULL, 
                          theme_set="theme_bw", theme_basesize=15,
                          x.ticks=NULL, x.ticks.lab=NULL,
                          x.ticks.angle=0, x.ticks.kb=TRUE,
                          x.gridmajor=TRUE, center.name="center", 
                          ylab="normalized counts", title=waiver(), 
                          subtitle=waiver(), colors=NULL, loess.span = NULL,
                          legend.ncol=NULL, lwd=.75, return.melted=FALSE,
                          CI.level=NULL, CI.reps=500,
                          CI.workers=NULL, CI.alpha=.2)

{
  
  #/ Checks:
  invisible(match.arg(class(names), c("character", "NULL")))
  invisible(match.arg(class(theme_set), c("character", "NULL")))
  invisible(match.arg(class(theme_basesize), c("numeric", "NULL")))
  invisible(match.arg(class(x.ticks), c("numeric", "NULL")))
  invisible(match.arg(class(x.ticks.lab), c("numeric", "NULL")))
  invisible(match.arg(class(x.ticks.angle), c("numeric")))
  invisible(match.arg(class(x.ticks.kb), c("logical")))
  invisible(match.arg(class(loess.span), c("numeric", "NULL")))
  invisible(match.arg(class(lwd), "numeric"))
  invisible(match.arg(class(CI.level), c("numeric", "NULL")))
  invisible(match.arg(class(CI.reps), c("numeric", "NULL")))
  invisible(match.arg(class(CI.workers), c("numeric", "NULL")))
  invisible(match.arg(class(CI.alpha), c("numeric", "NULL")))
  invisible(match.arg(class(return.melted), c("logical")))
  
  if(!is.null(CI.level)) if(!requireNamespace("boot",quietly=TRUE))
    stop("For CIs please install the boot package from CRAN", call. = FALSE)
  
  #/ check all ranges identical width
  if(class(ranges)=="GRanges") {
    w<-width(ranges)
  } else w<- unlist(width(ranges))
  if(var(w)!=0) stop("All intervals in ranges must have equal width") else w <- w[1]
  
  #/ check x.ticks not outside interval range
  if(!is.null(x.ticks)){
    if(sum(x.ticks < 1 | x.ticks > w)>0) {
      message("[Error] x.ticks are out of bounds. Allowed ticks range from 1 to ", w)
      message("        Offending ticks are ", 
              paste(x.ticks[x.ticks < 1 | x.ticks > w], collapse = ","))
      stop_quietly()
    }}
  
  #/ Set of colors from dittoSeq::dittoColors()
  cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
            "#0072B2", "#D55E00", "#CC79A7", "#666666", "#AD7700", 
            "#1C91D4", "#007756", "#D5C711", "#005685", "#A04700", 
            "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", 
            "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", 
            "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", 
            "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", 
            "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666", 
            "#3D3D3D")
  
  #############################
  # Create ScoreMatrix
  #############################
  
  #/ If multiple bigwigs
  if(length(bigwigs>1) & class(ranges)[1] %in% c("GRanges")){
    if(is.null(names)) names <- gsub(".bigwig|.BigWig|.bw|.BW", "",basename(bigwigs))
    if(length(bigwigs) != length(names)) stop("names is not of same length as bigwigs", call. = FALSE)
    sm.list <- lapply(1:length(bigwigs), function(x){
      genomation::ScoreMatrix(target=bigwigs[x], windows=ranges)
    })
  }
  
  #/ If one bigwig and one GRanges or a GRangesList
  if(length(bigwigs)==1 & class(ranges)[1] %in% "CompressedGRangesList"){
    if(is.null(names)){
      if(is.null(names(ranges))){
        names <- paste0("set", 1:length(ranges))
      } else names <- names(ranges)
    } else {
        if(length(ranges) != length(names)) 
          stop("names is not of same length as length of the GRangesList", call. = FALSE)
    }
    sm.list <- lapply(1:length(ranges), function(x){
      if(class(ranges)[1] %in% "GRanges") {
        rn <- ranges
      } else rn <- ranges[[x]]
      genomation::ScoreMatrix(target=bigwigs, windows=rn)
    })
  }
  
  if(length(bigwigs)>1 & class(ranges)[1] %in% "CompressedGRangesList"){
    stop("Use either one bigwig with a GRanges(list) or many bigwigs with a single GRanges")
  }
  
  #############################
  # Get colMeans
  #############################
  
  df.colmeans <- lapply(sm.list, function(x){
    data.frame(colMeans = base::colMeans(x))
  }) %>% do.call(cbind, .) %>% data.frame(Seq=1:nrow(.), .)
  colnames(df.colmeans)[2:ncol(df.colmeans)] <- names      
  
  #############################
  # Loess (optional)
  #############################
  
  if(!is.null(loess.span)){
    
    df.colmeans[,2:ncol(df.colmeans)] <- lapply(2:ncol(df.colmeans), function(x){
      
      d <- df.colmeans[,c(1,x)]; colnames(d)<-c("a","b")
      stats::loess(formula=b~a, data = d, span = loess.span)$fitted
      
    }) %>% do.call(cbind, .)
    
  }
  
  #############################
  # Assemble plot body / colors
  #############################
  
  df.melt <- reshape2::melt(df.colmeans, value.name = "colmeans",
                            id.vars = "Seq", variable.name = "samples")
  
  if(return.melted) return(df.melt)
  
  if(!is.null(theme_set)) theme_set(get(theme_set)(base_size=theme_basesize))
  
  #/ colors, either the user-provided or in-built ones:
  if(is.null(colors)){
    colors <- cols[1:length(sm.list)]
  } else {
    if(length(colors) != length(sm.list)){
      message("colors is not same length as curves to be plotted, using default colors")
      colors <- cols[1:length(sm.list)]
    }
  }
  
  p <- ggplot(df.melt, aes(x=Seq, y=colmeans, color=samples)) + 
    geom_line(size=lwd) +
    scale_color_manual(values=colors)
  
  #########################
  # Confidence intervals (optional)
  #########################
  
  if(!is.null(CI.level)){
    
    require(boot)
    
    if(is.null(CI.workers)) CI.workers <- detectCores()/2
    
    #/ Estimate confidence of the colMeans via CIs/boot:
    ci1 <- lapply(1:length(sm.list), function(x) 
      t(CalculateBootCI(cts = sm.list[[x]], Times = CI.reps, 
                        Cores = CI.workers, CI = CI.level))) %>%
      do.call(rbind, .)
    
    #/ Add to plot as ribbon:
    p <- p + geom_ribbon(aes(x = Seq, ymin = ci1[,1], ymax = ci1[,2], fill=samples), 
                         alpha=CI.alpha, show.legend = FALSE, colour = NA) +
      scale_fill_manual(values=colors)
    
  }
  
  #########################
  # Ticks and axis
  #########################
  
  #/ Default is to simple regard center as 0bp and then half the width 
  #/ to up/downstream, so a 2kb interval would be -1kb,0,+1kb
  if(is.null(x.ticks)){
    
    x.ticks <- c(0, ceiling(nrow(df.colmeans)/2), nrow(df.colmeans))
    if(x.ticks.kb) {
      bpkb<-"kb" 
      x.ticks.lab <- unlist(lapply(x.ticks,function(x) x/1000))
    } else {
        bpkb<-"bp" 
        x.ticks.lab<-x.ticks
    }
      
    x.ticks.lab <- paste0(x.ticks.lab,bpkb)
    x.ticks.lab[1] <- paste0("-",x.ticks.lab[2])
    x.ticks.lab[3] <- paste0("+",x.ticks.lab[2])
    x.ticks.lab[2] <- center.name
        
  }
  
  if(!is.null(x.ticks) & is.null(x.ticks.lab)) x.ticks.lab <- x.ticks
  
  #/ default legend behaviour, do rowwise filling with 3 elements per column
  if(is.null(legend.ncol)){
    legend.ncol=3
    if(length(names)==4) legend.ncol<-2
    legend.row=length(names)/2
  } else legend.row=NULL
  
  #/ add to plot, add legend:
  p <- p +
    labs(x = "", y = ylab) +
    scale_x_continuous(breaks=x.ticks,
                       labels=x.ticks.lab) +
    guides(x =  guide_axis(angle = x.ticks.angle)) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.key.width = unit(1.5,"cm"),
          legend.spacing.y = unit(0.5, 'cm')) +
    guides(col = guide_legend(ncol=legend.ncol,nrow=legend.row,byrow=TRUE)) +
    ggtitle(label = title, subtitle = subtitle)
  
  if(!x.gridmajor) 
    p <- p + theme(panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_blank()) 
  
  return(p)
  
}


