#' Perform genome segmentation depending
#' 
#' @param data A \code{GRanges} object containing the IP reads. See details for
#'  more information on how to set up the data.
#' @param dataC A \code{GRanges} object containing the control reads. Set to NULL
#'  by default, i.e. no control.
#' @param map A \code{GRanges} object containing the mappability profiles. Set
#'  to NULL by default, i.e. no profiles.
#' @param minReads A \code{numeric}. The minimum number of F/R reads to be 
#' present in the sliding window.
#' @param minReadsInRegion A \code{numeric}. The minimum number of F/R reads
#'  to be present in the region.
#' @param jitter A \code{logical} value stating whether some noise should be 
#' added to the read locations. This is recommended if the read positions have
#'  lots of duplicates.
#' @param maxLregion A \code{numeric}. The maximum length.
#' @param minLregion A \code{numeric}. The minimum length.
#' @param step A \code{numeric}. The increment of the sliding window.
#' @param width A \code{numeric}. The width of the region.
#' @param package A \code{character}. "PICS" or "PING"
#' 
#' @importFrom GenomicRanges seqnames strand
#' @export
segReadsGeneric <- function(data, dataC = NULL, map = NULL, minReads = 2, minReadsInRegion = 3,
                            jitter = FALSE, maxLregion = 0, minLregion = 100, step = 20,
                            width = 250, package = "PICS") {
    ## Check that we have the right data type
    if (!is(data, "GRanges")) {
        stop("The input data should be 'GRanges' object. Provided: ", class(data))
    }
    ## Check that we have the same number of chromosomes
    if (!is.null(dataC)) {
        if (length(levels(seqnames(dataC))) != length(levels(seqnames(data)))) {
            stop("Your IP and control data do not have the same number of chromosomes. IP: ", length(levels(seqnames(dataC))), " Control: ", length(levels(seqnames(data))))
        }
    }
    
    ## Total number of reads per sample
    lIP <- length(data)
    
    if (is.null(minReads)) {
        # Done once per chr
        chrs <- levels(seqnames(data))
        chrlist <- vector("list", length(chrs))
        names(chrlist) <- chrs
        for (cc in chrs) {
            chrlist[[cc]] <- diff(c(min(start(data[strand(data) == "+"])[1], end(data[strand(data) == "-"])[1]), max(tail(start(data[strand(data) == "+"]), 1), tail(end(data[strand(data) == 
                "-"]), 1))))
        }
        minReads <- (ceiling(as.numeric(lIP)/(2 * as.numeric(chrlist)) * width))
        minReads <- as.integer(names(which.max(table(minReads))))
        minReads <- min(minReads, 5)
        minReads <- max(minReads, 2)
        print(paste("We automatically calculated minReads, which is ", minReads, ".", sep = ""))
    }
    
    minReadsInRegion = max(minReadsInRegion, minReads)
    
    paraSW <- list(step = as.integer(step), width = as.integer(width), minReads = as.integer(minReads))
    if (!is.null(map) & !is(map, "GRanges")) {
        stop("Map should be a 'GRanges' object. Provided:", class(map))
    } else if (is.null(map)) {
        start <- NULL
        end <- NULL
    } else {
        map <- map[as.character(seqnames(map)) %in% levels(seqnames(data))]
        chrs <- levels(seqnames(data))
        start <- end <- vector("list", length(chrs))
        names(start) <- names(end) <- chrs
        for (cc in chrs) {
            start[[cc]] <- start(map[seqnames(map) == cc])
            end[[cc]] <- end(map[seqnames(map) == cc])
        }
    }
    
    if (maxLregion > 0) 
        maxStep = (maxLregion - 2 * paraSW$width)/paraSW$step else maxStep = 0
    
    ## Prepare C input:
    data <- .formatCInput(data)
    if (!is.null(dataC)) {
        lCont <- length(dataC)  #before data transformation
        dataC <- .formatCInput(dataC)
    } else {
        dataC <- vector("list", length(data))
        names(dataC) <- names(data)
        for (cc in names(data)) {
            dataC[[cc]] <- vector("list", 2)
            names(data[[cc]]) <- c("+", "-")
        }
        lCont <- 0
    }
    
    ## Perform the segmentation
    newSegReadsList <- .Call("segReadsAll", data, dataC, start, end, as.integer(jitter), paraSW, as.integer(maxStep), as.integer(minLregion), pPackage = package, 
        PACKAGE = "PICS")
    
    newSegReadsList <- unlist(newSegReadsList, recursive = FALSE, use.names = FALSE)
    if (is.null(newSegReadsList)) {
        stop("No Candidate regions found, you should decrease 'minReads'")
    }
    # newSet<-segReadsList(newSegReadsList,paraSW,as.integer(sum(unlist(lIP))),as.integer(sum(unlist(lCont))))
    newSet <- segReadsList(newSegReadsList, paraSW, as.integer(lIP), as.integer(lCont))
    
    ttt = summarySeg(newSet)
    indrm = ((ttt$L < minLregion) | (ttt$NF < minReadsInRegion) | (ttt$NR < minReadsInRegion))
    newSet@List = newSet@List[!indrm]
    
    return(newSet)
}

#' Summarize segmentList objects
#' 
#' Summarize segmentList objects into a data.frame
#' 
#' @param seg A \code{segmentList} object as returned by \code{segmentPICS}.
#' 
#' @return A \code{data.frame}. With 
#'  \itemize{
#'    \item{chr:} Chromosome id
#'    \item{NF:} Number of forward reads
#'    \item{NR:} Number of reverse reads
#'    \item{L:} Length of segment
#'    \item{min:} Start location of segments
#'    \item{max:} End location of segments
#'  }
#'  
#' @export
summarySeg <- function(seg) {
    temp <- .Call("getSegL", seg@List, PACKAGE = "PICS")
    ans <- data.frame(chr = temp[[1]], NF = temp[[2]], NR = temp[[3]], L = temp[[4]], min = temp[[5]], max = temp[[6]])
    ans$chr <- as.character(ans$chr)
    return(ans)
}


## Input: a GRanges object Output: a list: list$chr$strand
#' @importFrom IRanges start end
.formatCInput <- function(GRObject) {
    chrs <- levels(seqnames(GRObject))
    lData <- vector("list", length(chrs))
    names(lData) <- chrs
    for (cc in chrs) {
        GRccObject <- GRObject[seqnames(GRObject) == cc]
        lData[[cc]] <- vector("list", 2)
        names(lData[[cc]]) <- c("+", "-")
        lData[[cc]][["+"]] <- start(GRccObject[strand(GRccObject) == "+"])
        lData[[cc]][["-"]] <- end(GRccObject[strand(GRccObject) == "-"])
    }
    return(lData)
}
