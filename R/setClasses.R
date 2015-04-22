#' segReads class
#' 
#' This class stores a candidate region, as determined by a segmentation function.
#' It contains information regarding the reads that contributes to this region.
#' Usually, \code{segReads} are part of a \code{segReadsList} object, to be used
#' as input for the \code{PICS} and \code{PING} function.
#' 
#' @slot yF A \code{numeric}. The start position of forward reads in the 
#'  candidate region.
#' @slot yR A \code{numeric}. The start position of reverse reads in the 
#'  candidate region.
#' @slot cF A \code{numeric}. The start position of forward reads in the control.
#' @slot cR A \code{numeric}. The start position of reverse reads in the control.
#' @slot map A \code{matrix}. Mappability profile
#' @slot chr A \code{character}. The chromosome on which the candidate region is
#'  located.
#'  
#' @author Xuekui Zhang, Arnaud Droit, Raphael Gottardo, Renan Sauteraud
#' @seealso \code{\linkS4class{segReadsList}}, \code{\link{segmentPICS}}, \code{\link{PICS}}
#' @aliases
#' segReads
#' show,segReads-method
#' @export
setClass("segReads",
         representation(yF = "numeric", yR = "numeric", cF = "numeric",
                        cR = "numeric", map = "matrix", chr = "character"),
         prototype(yF = numeric(0), yR = numeric(0), cF = numeric(0),
                   cR = numeric(0), map = matrix(0, 0, 2), chr = character(0)))

segReads<-function(yF, yR, cF, cR, map, chr){
  if(!is.vector(yF) || !is.vector(yR) || !is.numeric(yF) || !is.numeric(yR)){
    stop("Argument 'yF/yR' must be numeric vectors ", call. = FALSE)
  }
  if((!is.vector(cF) || !is.vector(cR) || !is.numeric(cF) || !is.numeric(cR)) & (!is.null(cF) || !is.null(cR))){
    stop("Argument 'cF/cR' must be numeric vectors ", call. = FALSE)
  }
  if(!is.matrix(map)){
    stop("Argument 'map' must be a matrix ", call. = FALSE)
  }
  if(!is.character(chr)){
    stop("Argument 'chr' must be a character string", call. = FALSE)
  }
  new("segReads", yF=yF, yR=yR, cF=cF, cR=cR, map=map, chr=chr)
}


#' segReadsList class
#' 
#' This class stores the results of \code{segmentPICS}, the list of all candidate
#' regions to be used as input in \code{PICS} or \code{PING}. On top of
#' holding all the \code{segReads}, it contains additional information regarding
#' the parameters used to generate this object.
#' 
#' @slot List A \code{list} of \code{segReads} objects. The list of candidate
#'  regions.
#' @slot paraSW A \code{list} of \code{numeric}. A list of parameters used in
#'  the segmentation function.
#' @slot N A \code{numeric}. The number of reads used as input of the
#'  segmentation function.
#' @slot Nc A \code{numeric}. The number of control reads used as input of the
#'  segmentation function.
#' 
#' @note
#' The segmentation functions are \code{segmentPICS} and \code{segmentPING}
#' 
#' @section Methods:
#' \describe{
#'  \item{show}{signature(object = "segReadsList"): Print details about the
#'   object.}
#'  \item{length}{signature(x = "segReadsList"): The length of the \code{List}
#'   in the object. The number of candidate regions identified.}
#' } 
#' 
#' 
#' @author Xuekui Zhang, Arnaud Droit, Raphael Gottardo, Renan Sauteraud
#' @seealso \code{\linkS4class{segReads}}, \code{\link{segmentPICS}}, \code{\link{PICS}}
#' @aliases
#' segReadsList
#' show,segReadsList-method
#' length,segReadsList-method
#' [,segReadsList,ANY,ANY-method
#' [[,segReadsList,ANY,ANY-method
#' @export
setClass("segReadsList",
         representation(List = "list", paraSW = "list",
                        N = "integer", Nc = "integer"),
         prototype(List = list(List = list(0),
                               paraSW = list(step = integer(0), width = integer(0), minReads = integer(0)),
                               N = integer(0), Nc = integer(0))))

segReadsList<-function(List,paraSW,N,Nc){
  if(!is.list(paraSW) & !all(sapply(paraSW,"is.numeric"))){
    stop("Argument 'paraSW' must be a list of numeric arguments", call.=FALSE)
  }
  if(any(lapply(List,"class")!="segReads")){
    stop("Argument 'List' must be a list of segReads arguments", call.=FALSE)
  }
  if(!is.integer(N) | !is.integer(Nc)){
    stop("Argument 'N' and 'Nc' must be integers", call.=FALSE)    
  }
  new("segReadsList", List=List, paraSW=paraSW, N=N, Nc=Nc)
}

#' pics class
#' 
#' This class store the result of a \code{PICS} analysis: The estimated location
#' of a binding event, its score and the chromosome it belongs to.
#' 
#' @slot estimates A \code{list} of \code{numeric}. Coordinates of the binding
#'  event, variance and standard error. Contains the following elements:
#' \describe{
#'   \item{w}{The mixture weights. The proportion of reads of the region 
#'   contributing to each binding event.}
#'   \item{mu}{The binding site positions.}
#'   \item{delta}{The average DNA fragment lengths. It corresponds to the
#'    distance between the maxima of the forward and reverse read position
#'    densities.}
#'   \item{sigmaSqF}{The variance parameters for the forward distribution.}
#'   \item{sigmaSqR}{The variance parameters for the reverse distribution.}
#'   \item{seMu}{The standard errors for mu.}
#'   \item{seMuF}{The standard errors for muF.}
#'   \item{seMuR}{The standard errors for muR.}
#' }
#' @slot score A \code{numeric}. The score for the binding event.
#' @slot scoreF A \code{numeric}. The score of the forward reads for the binding
#'  event.
#' @slot scoreR A \code{numeric}. The score of the reverse reads for the binding
#'  event.
#' @slot Nmerged A \code{numeric}. The number of peaks that got merged (in case
#'  of boundary issue).
#' @slot converge A \code{logical}. TRUE if the EM has converged.
#' @slot range A \code{numeric}. The range of the binding event.
#' @slot chr A \code{character}. The chromosome where the binding event is located.
#'  
#' @author Xuekui Zhang, Arnaud Droit, Raphael Gottardo, Renan Sauteraud
#' @seealso \code{\linkS4class{picsList}}, \code{\link{PICS}}
#' @aliases
#' show,pics-method
#' w
#' w,pics-method
#' mu
#' mu,pics-method
#' delta
#' delta,pics-method
#' se
#' se,pics-method
#' seF
#' seF,pics-method
#' seR
#' seR,pics-method
#' sigmaSqF
#' sigmaSqF,pics-method
#' sigmaSqR
#' sigmaSqR,pics-method
#' score,pics-method
#' scoreForward
#' scoreForward,pics-method
#' scoreReverse
#' scoreReverse,pics-method
#' chromosome
#' chromosome,pics-method
#' as.data.frame.picsList
#' @export
setClass("pics",
         representation(estimates = "list", score = "numeric", scoreF = "numeric",
                        scoreR = "numeric", Nmerged = "numeric",
                        converge ="logical", range = "numeric", chr = "character"),
         prototype(estimates = list(w = numeric(0), mu = numeric(0), 
                                    delta = numeric(0), sigmaSqF = numeric(0),
                                    sigmaSqR = numeric(0), seMu = numeric(0),
                                    seMuF = numeric(0), seMuR = numeric(0)),
                  score = numeric(0), scoreF = numeric(0), scoreR = numeric(0),
                  Nmerged = numeric(0), converge = logical(0),
                  range = numeric(0), chr = character(0)))


#' picsError class
#' 
#' This object is used to return an error code when the PICS function failed to
#' return a valid set of estimates for a candidate regions. This could be due to
#' non-convergence of the EM algorithm, a singular information matrix, or a
#' number of reads below the limit specified by the user. All of these are
#' typically due to too few reads in the region and do not affect the rest of
#' the analysis, as such regions would most likely be labelled as false
#' positives. 
#' 
#' @author Xuekui Zhang, Arnaud Droit, Raphael Gottardo, Renan Sauteraud
#' @seealso \code{\linkS4class{picsList}}, \code{\link{PICS}}
#' @aliases
#' show,picsError-method
#' score,picsError-method
#' scoreForward,picsError-method
#' scoreReverse,picsError-method
#' chromosome,picsError-method
#' w,picsError-method
#' mu,picsError-method
#' delta,picsError-method
#' sigmaSqF,picsError-method
#' sigmaSqR,picsError-method
#' se,picsError-method
#' seF,picsError-method
#' seR,picsError-method
#' @export
setClass("picsError",
         representation(errorCode = "character"),
         prototype(errorCode = character(0)))


#' picsList class
#' 
#' This class stores the result of a \code{PICS} analysis. It consists in a list
#' of all binding events identified held by \code{pics} objects and some 
#' additional information regarding the parameters used to generate this object.
#' 
#' @slot List A \code{list} of \code{pics} objects. The list of predicted
#'  binding events.
#' @slot paraPrior A \code{list} of \code{numeric}. The list of prior parameters.
#' @slot paraEM A \code{list}. The list of parameters used in the EM algorithm
#'  of the \code{PICS} analysis.
#' @slot minReads A \code{list} of \code{integer}. The minumum number of reads
#'  per peak and per region.
#' @slot N A \code{numeric}. The number of reads used in the analysis.
#' @slot Nc A \code{numeric}. The number of control reads used in the analysis.
#'
#' @section Methods:
#' \describe{
#'  \item{show}{signature(object = "segReadsList"): Print details about the
#'   object.}
#'  \item{length}{signature(x = "segReadsList"): The length of the \code{List}
#'   in the object. The number of candidate regions identified.}
#' }
#'  
#' @author Xuekui Zhang, Arnaud Droit, Raphael Gottardo, Renan Sauteraud
#' @seealso \code{\linkS4class{pics}}, \code{\linkS4class{picsError}}, \code{\link{PICS}}
#' @aliases
#' picsList
#' show,picsList-method
#' length,picsList-method
#' plot,picsList,picsList-method
#' score,picsList-method
#' scoreForward,picsList-method
#' scoreReverse,picsList-method
#' chromosome,picsList-method
#' w,picsList-method
#' mu,picsList-method
#' delta,picsList-method
#' sigmaSqF,picsList-method
#' sigmaSqR,picsList-method
#' se,picsList-method
#' seF,picsList-method
#' seR,picsList-method
#' [,picsList,ANY,ANY-method
#' [[,picsList,ANY,ANY-method
#' @export
setClass("picsList",
         representation(List = "list", paraPrior = "list", paraEM = "list",
                        minReads = "list", N = "integer", Nc = "integer"),
         prototype(List = list(0), minReads = list(perPeak = integer(0), perRegion = integer(0)),
                   paraPrior = list(xi = double(0), rho = double(0),
                                    alpha = double(0), beta = double(0)),
                   paraEM = list(kMax = integer(0), B = integer(0), tol = double(0)),
                   N = integer(0), Nc = integer(0)))

newPicsList <- function(List, paraEM, paraPrior, minReads, N, Nc){
  if(!is.list(paraEM) & !all(sapply(paraEM,"is.numeric"))){
    stop("Argument 'paraEM' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(paraPrior) & !all(sapply(paraPrior,"is.numeric"))){
    stop("Argument 'paraPrior' must be a list of numeric arguments", call.=FALSE)
  }
  if(!is.list(minReads) & !all(sapply(minReads,"is.numeric"))){
    stop("Argument 'minReads' must be a list of numeric arguments", call.=FALSE)
  }
  if(!all((lapply(List,"class")=="pics" | lapply(List,"class")=="picsError"))){
    stop("Argument 'List' must be a list of 'pics' or 'picsError' arguments", call.=FALSE)
  }
  if(!is.integer(N) | !is.integer(Nc)){
    stop("Argument 'N' and 'Nc' must be integers", call.=FALSE)    
  }
  new("picsList", List=List, paraEM=paraEM, paraPrior=paraPrior,
      minReads=minReads, N=N, Nc=Nc)
}