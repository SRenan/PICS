#' makeGRanges
#' 
#' Generate a \code{GRanges} objects from a \code{picsList} object.
#' 
#' @param obj A \code{picsList} object, as returned by \code{PICS}.
#' @param type A \code{character}. The type of intervals to be created. The
#'  valid types are  "bed', "wig" and "fixed". See details section for more.
#' @param filter A \code{list} of filters to be used before computing the FDR.
#'  By default, all regions are included.
#' @param length A \code{numeric}. The length of the intervals when using 
#'  "fixed" type.
#'
#' @details
#' "bed" will generate intervals from the forward peak max to the reverse peak 
#' max. "wig" will generate a density profile for the forward and reverse reads.
#' "bed" and "wig" types should be used to be exported to wig/bed files that can
#' be used with the UCSC genome browser. "fixed" corresponds to the binding site
#' estimates +/-3*length. "bed" and "wig" files can be exported using the
#' \code{export} function of the \code{rtracklayer} package.
#' 
#' 
#' @return A \code{GRanges} object
#' 
#' @author Renan Sauteraud
#' @seealso \code{\link{picsFDR}}
#' 
#' @importClassesFrom IRanges IRanges
#' @importFrom IRanges IRanges
#' @importClassesFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRanges findOverlaps
#' @export
makeGRanges <- function(obj, type = "fixed",
                        filter = list(delta = c(0,Inf), se = c(0,Inf), 
                                      sigmaSqF = c(0,Inf), sigmaSqR = c(0,Inf),
                                      score = c(0,Inf)), length = 100){
  nSe<-3
  if(type!="wig"){
    mu      <- mu(obj)
    delta   <- delta(obj)
    se      <- se(obj)
    seF     <- seF(obj)
    seR     <- seR(obj)
    sf      <- sigmaSqF(obj)
    sr      <- sigmaSqR(obj)
    chromosome <- chromosome(obj)
    score   <- score(obj)
    
    ## Filter regions with small deltas over a region
    if(!is.null(filter)){
      # indK<-unlist(sapply(obj@List,function(x,dmin){k<-K(x);if(k==0){return(NULL);};if(any(delta(x)<dmin & se(x)<20)){return(rep(FALSE,k));}else{return(rep(TRUE,k));}},dmin=filter$delta[1]))
      ### Filter based on delta
      ind1<-delta>filter$delta[1] & delta<filter$delta[2]
      ind2<-sf>filter$sigmaSqF[1] & sf<filter$sigmaSqF[2]
      ind3<-sr>filter$sigmaSqR[1] & sr<filter$sigmaSqR[2]
      ind5<-is.finite(score) & score>filter$score[1] & score<filter$score[2]
	    ind<-ind1&ind2&ind3&ind5
      
      if(!is.null(filter$se)){
      	ind4<-se>filter$se[1] & se<filter$se[2] & seF>filter$se[1] & seF<filter$se[2] & seR>filter$se[1] & seR<filter$se[2]
      	ind<-ind&ind4
      }
    } else{
      ind<-is.finite(score)
    }
  }
  if(type=="bed"){
    score<-(score(obj))[ind]
    ord<-order(-score)
    # Order the score    
    score<-score[ord]
    start<-(mu-delta/2-nSe*seF)[ind]
    end<-(mu+delta/2+nSe*seR)[ind]
    start<-start[ord]
    end<-end[ord]
    # chrom<-(paste("chr", chromosome, sep=""))[ord]    
    chrom<-(paste("",chromosome, sep=""))[ind]
    chrom<-chrom[ord]
    # strand<-NULL
  } else if(type=="ci"){
    score<-(score(obj))[ind]
    ord<-order(-score)
    score<-score[ord]
    start<-(mu-nSe*se)[ind]
    end<-(mu+nSe*se)[ind]
    start<-start[ord]
    end<-end[ord]
    chrom<-(paste("",chromosome, sep=""))[ind]
    chrom<-chrom[ord]
    # strand<-NULL
  } else if(type=="fixed"){
    score<-(score(obj))[ind]
    ord<-order(-score)
    score<-score[ord]
    start<-(mu-length)[ind]
    end<-(mu+length)[ind]
    start<-start[ord]
    end<-end[ord]
    chrom<-(paste("",chromosome, sep=""))[ind]
    chrom<-chrom[ord]
    strand<-NULL
  } else if(type=="wig"){
    temp<-wigDensity(obj,strand="*",step=10,sum=TRUE,filter=filter,scale=TRUE)
    chrom<-temp$chr
    start<-temp$x
    end<-temp$x+9
    score<-temp$density
    strand="*"
  }
  
  if(length(score)==0) {
  	gr=NULL
  } else{
	  ranges<-IRanges(as.integer(start),as.integer(end))
	  if(type=="bed" | type=="ci" | type=="fixed"){
		  names(ranges)<-paste("pics",1:(length(score)),sep="")
		  gr<-GRanges(seqnames = chrom, ranges = ranges, strand = "*", score)
	  } else{ #wig
		  gr<-GRanges(seqnames = chrom, ranges = ranges, strand=strand, score)
		  
		  print("Removing overlapping binding events")
		  overlap<-as.matrix(findOverlaps(ranges(gr)))
		  overlap<-split(overlap[,1], overlap[,2])
		  scL<-score(gr)
		  maxScores <- lapply(overlap, function(x){
        max(scL[x])
        })
		  maxScscores<-as.numeric(unlist(maxScores))
		  idx<-which(maxScores==scL)
		  gr<-gr[idx,]
		  
#		  gr<-killOverlaps(gr)
		  #Create the new GRanges objec to return
#		  gr<-GRanges(seqnames=chrom, ranges=ranges, score, strand=strand)
	  }
  } 	
  return(gr)
}

wigDensity <- function(x, strand="+", step=10, sum=FALSE, filter=NULL, scale=TRUE){
  # Check that all filters are passed
  missingNames <- !c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score") %in% names(filter)
  filter[c("delta","sigmaSqF","sigmaSqR","se","seF","seR","score")[missingNames]]<-list(c(0,Inf))
  if(strand=="+"){
    strand <- 1
  }
  else if(strand=="-"){
    strand <- -1
  }
  else if(strand=="*"){
    strand<-0
  }
  else{
    stop("Strand must be either '+', '-' or '*'")
  }
  ans<-.Call("getDensityList", x, strand, step, filter, sum, scale, PACKAGE="PICS")
  return(ans)
}