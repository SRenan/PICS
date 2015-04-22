#' Segment the genome into candidate regions
#' 
#' Pre-process bidirectional aligned reads data from a single ChIP-Seq
#' experiment to detect candidate regions with a minimum number of forward and
#' reverse reads. These candidate regions will then be processed by PICS.
#' 
#' @aliases segmentPICS segReadsGeneric
#' @param data A file or object containing the IP reads.
#' See details for more information on how to set up the data.
#' @param dataC A file or object containing the control
#' reads. Set to NULL by default, i.e. no control.
#' @param map A file or object containing the mappability profiles. Set to
#' NULL by default, i.e. no profiles.
#' @param minReads The minimum number of F/R reads to be present in the sliding
#' window.
#' @param minReadsInRegion The minimum number of F/R reads to be present in the
#' region.
#' @param jitter A logical value stating whether some noise should be added to
#' the read locations. This is recommended if the read positions have lots of
#' duplicates.
#' @param maxLregion The maximum length.
#' @param minLregion The minimum length.
#' @return An object of class \code{segReadsList} containing the results for
#' all regions pre-processed.
#' 
#' @details
#' The input of data, dataC and map can be an object of class \code{GRanges} or 
#' \code{GAlignments}. Or a bed or bam file.
#' 
#' @author Renan Sauteraud
#' @seealso \code{\linkS4class{segReadsList}}
#' @references X. Zhang, G. Robertson, M. Krzywinski, K. Ning, A. Droit, S.
#' Jones, and R. Gottardo, PICS: Probabilistic Inference for ChIP-seq
#' arXiv, 0903.3206, 2009.
#' @keywords functions
#' @examples
#' 
#' # Read data
#' datadir <- system.file("extdata", package="PICS")
#' dataIP <- file.path(datadir, "Treatment_tags_chr21_sort.bed")
#' dataCont <- file.path(datadir, "Input_tags_chr21_sort.bed")
#' map <- file.path(datadir, "mapProfileShort.bed")
#' 
#' seg<-segmentPICS(dataIP, dataC = dataCont, map = map, minReads = 1)
#' 
#' @importFrom tools file_ext
#' @importFrom GenomicAlignments readGAlignmentsFromBam
#' @importClassesFrom GenomicAlignments GAlignments
#' @export
segmentPICS <- function(data, dataC = NULL, map = NULL, minReads = 2,
                      minReadsInRegion = 3, jitter = FALSE,
                      maxLregion = 0, minLregion = 100){
  
  step <- 20
  dataType <- "TF"
  if(dataType == "TF")  width <- 250
  if(dataType == "H")   width <- 150 #Histones not implemented yet
  
  data  <- .check_segmentPICS_input(data, "data")
  dataC <- .check_segmentPICS_input(dataC, "control")
  map   <- .check_segmentPICS_input(map, "map")
  
  newSet <- segReadsGeneric(data, dataC = dataC, map = map, minReads = minReads,
                            minReadsInRegion = minReadsInRegion, jitter = jitter,
                            maxLregion = maxLregion, minLregion = minLregion,
                            step = step, width = width, package = "PICS")
  return(newSet)
}

.check_segmentPICS_input <- function(data, type = "data"){
  if(is(data, "GAlignments")){
    data <- as(data, "GRanges")
  } else if(is(data, "character")){
    ext <- file_ext(data)
    if(ext == "bed"){
      data <- as(read.table(data, header = TRUE), "GRanges")
    } else if(ext == "bam"){
      data <- as(readGAlignmentsFromBam(data), "GRanges")
    } else{
      stop("The given file should be in bam or bed format.")
    }
  } else if(is(data, "GRanges")){
  } else if(is.null(data)){
    if(type == "data"){
      stop("The data argument cannot be NULL.")
    }
  } else{
    stop(paste("The", type, "should be an object of class 'Granges' or 'GAlignments' or
         a filename."))
  }
  return(data)
}
