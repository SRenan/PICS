## show ------------------------------------------------------------------------
#' @export
setMethod("show", "segReads", function(object){
  cat("Object of class ", as.character(class(object)), "\n")
  cat("This object has the following slots: \n")
  cat(paste(names(getSlots(class(object))), collapse=", "), "\n")
})

#' @export
setMethod("show", "segReadsList", function(object){
  cat("Object of class",as.character(class(object)), "\n")
  cat("This object has the following slots: \n")
  cat(paste(names(getSlots(class(object))), collapse=", "), "\n")
  cat(paste("List is a list of", length(object),
            "'segReads' ojects, each of which has the following slots:\n"))
  cat("yR, yF, cR, cF, map, chr\n")
})

#' @export
setMethod("show", "pics", function(object){
  cat("Object of class ", class(object), "\n")
  cat("This object has the following slots: \n")
  cat(paste(names(getSlots(class(object))), collapse=", "), "\n")
})

#' @export
setMethod("show", "picsError", function(object){
  cat("Object of class ", class(object), "\n")
  cat("This object has the following slot: \n")
  cat(paste(names(getSlots(class(object))), collapse=", "), "\n")   
})

#' @export
setMethod("show", "picsList", function(object){
  cat("Object of class ", class(object), "\n")
  cat("This object has the following slots: \n")
  cat(paste(names(getSlots(class(object))), collapse=", "), "\n")
  cat(paste("List is a list of", length(object),
            "'pics' or 'picsError' ojects\n"))
})

## coerce ----------------------------------------------------------------------
setAs("picsList", "data.frame", function(from){
  List <- from@List
  List <- List[which(sapply(List, class) != "picsError")]
  ans <- data.frame(ID = rep(1:length(List), 
                             sapply(List, function(x){length(x@estimates$w)})),
                    w = w(from), mu = mu(from), delta = delta(from),
                    sigmaSqF = sigmaSqF(from), sigmaSqR = sigmaSqR(from),
                    se = se(from), score = score(from),
                    scoreF = scoreForward(from), scoreR = scoreReverse(from),
                    chr = chromosome(from))
  return(ans)                                                             
})
#' @export
as.data.frame.picsList <- function(x, row.names = NULL, optional =FALSE, ...){
  res <- as(x, "data.frame")
  return(res)
}

## length ----------------------------------------------------------------------
#' @export
setMethod("length", "segReadsList", function(x){
  return(length(x@List))
})

#' @export
setMethod("length", "picsList", function(x){
 return(length(x@List))
})

## subset ----------------------------------------------------------------------
#' @export
setMethod("[", "segReadsList", function(x,i, j,..., drop=FALSE){
	if(missing(i)){
		return(x)
	}
	if(!missing(j)){
		  stop("incorrect number of dimensions")
	} else{
    segReadsList(x@List[i], x@paraSW, x@N, x@Nc)
  }
})

#' @export
setMethod("[[","segReadsList", function(x, i, j, ..., exact = TRUE){
  if(length(i) != 1){
    stop("subscript out of bounds (index must have length 1)")
  }
  if(missing(i)){
    return(x)
  }
  if(!missing(j)){
    stop("incorrect number of dimensions")
  }
  x@List[[i]]
})

#' @export
setMethod("[","picsList",	function(x,i, j,..., drop=FALSE){
	if(missing(i)){
		return(x)
	}
	if(!missing(j)){
	  stop("incorrect number of dimensions")
	} else{
    newPicsList(x@List[i], x@paraEM, x@paraPrior, x@minReads, x@N, x@Nc)        
  }
})

#' @export
setMethod("[[","picsList", function(x, i, j, ..., exact = TRUE){
  if(length(i) != 1){
    stop("subscript out of bounds (index must have length 1)")
  }
  if(missing(i)){
    return(x)
  }
  if(!missing(j)){
    stop("incorrect number of dimensions")
  }
  x@List[[i]]
})

## plot ------------------------------------------------------------------------
#' @importFrom graphics grid axis mtext abline par
#' @export
setMethod("plot", signature("picsList", "picsList"),
          function(x, y, filter=NULL, h=.1, ...){
  FDR <- picsFDR(x, y, filter = filter)
  arg <- list(...)
  par(mar = c(4, 4, 4.5, 4) + 0.1)
  # points(FDR[,2],FDR[,3]/max(FDR[,3]),xaxt="n",yaxt="n",lty=3,col=3,pch=2)
  if(length(arg$xlim) != 2){
    xlim<-range(FDR[,2])
  }
  else{
    xlim <- c(max(arg$xlim[1], min(FDR[,2])), min(arg$xlim[2], max(FDR[,2])))
  }
  plot(FDR[,2], FDR[,1], xlab="score", ylab="FDR", panel.first=grid(nx=50), ...)
  xx <- FDR[FDR[,2]>xlim[1] & FDR[,2]<xlim[2],2]
  yy <- FDR[FDR[,2]>xlim[1] & FDR[,2]<xlim[2],3]
  xx <- xx[seq(1,length(xx), length.out=10)]
  yy <- yy[seq(1,length(yy), length.out=10)]
  axis(3, at=xx, labels=yy)
  mtext("# regions", side = 3, line = 3)
  FDRex <- FDR[FDR[,1]>0,]
  notDup <- rev(!duplicated(rev(FDRex[,1])))
  lines(FDRex[notDup,2], FDRex[notDup,1], col=2, lty=2, lwd=1.5)
  abline(h=h, lw=1.5, col="grey")  
})


## accessors -------------------------------------------------------------------
#' @export
setGeneric("w", function(x, ...) standardGeneric("w"))
#' @export
setMethod("w", "pics", function(x){
  return(x@estimates$w)
})
#' @export                                                                
setMethod("w", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("w", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(0), PACKAGE="PICS")
  return(ans)
})

#' @export
setGeneric("mu", function(x, ...) standardGeneric("mu"))
#' @export
setMethod("mu", "pics", function(x){
  return(x@estimates$mu)
})
#' @export                                                                
setMethod("mu", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("mu", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(1), PACKAGE="PICS")
  return(ans)
})

#' @export
setGeneric("delta", function(x, ...) standardGeneric("delta"))
#' @export
setMethod("delta", "pics", function(x){
  return(x@estimates$delta)
})
#' @export
setMethod("delta", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("delta", "picsList", function(x){
  ans <- .Call("getVector", x@List, as.integer(2), PACKAGE="PICS");
  return(ans)
})

#' @export
setGeneric("se", function(x, ...) standardGeneric("se"))
#' @export
setMethod("se", "pics", function(x){
  return(x@estimates$seMu)
})
#' @export
setMethod("se", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("se", "picsList", function(x){
  ans <- .Call("getVector", x@List, as.integer(5), PACKAGE="PICS");
  return(ans)
})

#' @export
setGeneric("sigmaSqF", function(x, ...) standardGeneric("sigmaSqF"))
#' @export
setMethod("sigmaSqF", "pics", function(x){
  return(x@estimates$sigmaSqF)
})
#' @export
setMethod("sigmaSqF", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("sigmaSqF", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(3), PACKAGE="PICS");
  return(ans)
})

#' @export
setGeneric("sigmaSqR", function(x, ...) standardGeneric("sigmaSqR"))
#' @export
setMethod("sigmaSqR", "pics", function(x){
  return(x@estimates$sigmaSqR)
})
#' @export
setMethod("sigmaSqR", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("sigmaSqR", "picsList", function(x){
  ans <- .Call("getVector", x@List, as.integer(4), PACKAGE="PICS");
  return(ans)
})

#' @importMethodsFrom IRanges score
#' @export
setMethod("score", "pics", function(x){
  return(x@score)
})
#' @importMethodsFrom IRanges score
#' @export
setMethod("score", "picsError", function(x){
  return(NULL)
})
#' @importMethodsFrom IRanges score
#' @export
setMethod("score", "picsList", function(x){
  ans <- .Call("getScore", x@List, PACKAGE="PICS")
  return(ans)
})
#' @export
setMethod("score", "picsError", function(x){
  return(NULL)                                                  
})
#' @export
setMethod("score", "picsList", function(x){
  ans<-.Call("getScore", x@List, PACKAGE="PICS")
  return(ans)
})

#' @export
setGeneric("scoreReverse", function(x, ...) standardGeneric("scoreReverse"))
setMethod("scoreReverse", "pics", function(x){
  return(x@scoreR)                                              
})                                                                        
#' @export
setMethod("scoreReverse", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("scoreReverse", "picsList", function(x){
  ans<-.Call("getScoreR", x@List, PACKAGE="PICS")
  return(ans)
})

#' @export
setGeneric("scoreForward", function(x, ...) standardGeneric("scoreForward"))
#' @export
setMethod("scoreForward", "pics", function(x){
  return(x@scoreF)                            
})
#' @export
setMethod("scoreForward", "picsError", function(x){
  return(NULL)
})
#' @export                                                                
setMethod("scoreForward", "picsList", function(x){
  ans<-.Call("getScoreF", x@List, PACKAGE="PICS")
  return(ans)
})

#' @export
setGeneric("chromosome", function(x, ...) standardGeneric("chromosome"))
#' @export
setMethod("chromosome", "pics", function(x){
  return(rep(x@chr, length(x@estimates$w)))
})
#' @export
setMethod("chromosome", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("chromosome", "picsList", function(x){
  ans<-.Call("getChr", x@List, PACKAGE="PICS")
  return(ans)
})

#' @export
setGeneric("seF", function(x, ...) standardGeneric("seF"))
#' @export
setMethod("seF", "pics", function(x){
  return(x@estimates$seMuF)
})
#' @export
setMethod("seF", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("seF", "picsList", function(x){
  ans <- .Call("getVector", x@List, as.integer(6), PACKAGE="PICS")
  return(ans)
})                                                                         

#' @export
setGeneric("seR", function(x, ...) standardGeneric("seR"))
#' @export
setMethod("seR", "pics", function(x){
  return(x@estimates$seMuR)
})
#' @export
setMethod("seR", "picsError", function(x){
  return(NULL)
})
#' @export
setMethod("seR", "picsList", function(x){
  ans<-.Call("getVector", x@List, as.integer(7), PACKAGE="PICS")
  return(ans)
})                                                                         
