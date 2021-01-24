#' Stop quietly
stop_quietly <- function() {
  # https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

#' Calculate column-wise confidence intervals from a numeric matrix
#'
#' Function to calculate confidence interval from a numeric matrix on a column-wise basis.
#' Code is adapted and basically identical to the Utilities.R script from the Bioconductor package ChIPSeeker.
#' Main intend is to serve as helper function for the PlotProfiles() function.
#' 
#' @param cts a count matrix
#' @param Times number of resamplings
#' @param Cores numbers of mclapply workers
#' @param CI the confidence interval, default is 95\%
#' 
#' @details 
#' 
#' @examples 
#' cts <- sapply(seq(1,100), function(x) rnorm(100, 1))
#' ci <- CalculateBootCI(cts = cts, Times = 500, Cores = 8, CI = 0.95)
#' 
#' @author Alexander Toenges
#' @export
CalculateBootCI <- function(cts, Times = 500, Cores = 1, CI = 0.95){
  
  if(!class(cts)[1] %in% c("matrix", "ScoreMatrix")) stop("cts must be a matrix or ScoreMatrix", call. = FALSE)
  if(Sys.info()['sysname'] %in% c("Darwin", "Linux")) {
    pl="multicore"
  } else pl="none"; Cores=1
  
  #/ The "statistic" function for creating the bootstrap replicates.
  #/ data is the original data and idx is the rows the original data are being subset to
  getcolsum <- function(data, idx) matrixStats::colMeans2(data[idx,])
  
  # parser:
  parseBootCiPerc <- function(bootCiPerc){
    bootCiPerc <- bootCiPerc$percent
    tmp <- length(bootCiPerc)
    ciLo <- bootCiPerc[tmp - 1]
    ciUp <- bootCiPerc[tmp]
    return(c(ciLo, ciUp))
  }
  
  #/ the main bootstrap function generating the replicates:
  booty <- boot::boot(data = cts, statistic = getcolsum, R = Times,
                      parallel = pl, ncpus = Cores)
  
  booty.CI <- 
    sapply(seq_len(ncol(cts)), 
           function(i) parseBootCiPerc(boot::boot.ci(booty, type = "perc", index = i, conf = CI)))
  
  return(booty.CI)
  
}