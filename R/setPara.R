## Functions to return a list of parameters to be used by PICS functions
#The default paraneters are set for PICS

#' Function that returns a list of parameters for the EM algorithm that can be
#' used as an argument of PICS.
#' 
#' This function takes from 0 to 7 EM algorithm parameters as argument, check
#' if they are valid and returns a list to be used in a call to PICS.
#' 
#' 
#' @param minK An \code{integer}. The minimum number of binding events per
#' region. If the value is 0, the minimum number is automatically calculated.
#' @param maxK An \code{integer}. The maximum number of binding events per
#' region. If the value is 0, the maximum number is automatically calculated.
#' @param tol A \code{numeric}. The tolerance for the EM algorithm.
#' @param B An \code{integer}. The maximum number of iterations to be used.
#' @param mSelect A character string specifying the information criteria to be
#' used when selecting the number of binding events.
#' @param mergePeaks A \code{logical} stating whether overlapping binding
#' events should be picked.
#' @param mapCorrect A \code{logical} stating whether mappability profiles
#' should be incorporated in the estimation, i.e: missing reads estimated.
#' @param dataType A \code{character}. If a dataType is set, the algorithm will
#' use the default parameters for this type of data (all the previous arguments
#' will be ignored). When using \code{PICS}, this argument should be NULL. For
#' \code{PING}, "MNase" and "Sonicated" are valid.
#' 
#' @return Returns a list of parameters to be used in PICS.
#' @author Renan Sauteraud
#' @seealso \code{PICS}
#' @keywords functions
#' @export setParaEM
setParaEM <- function(minK = 1, maxK = 15, tol = 1e-4, B = 100, mSelect = "BIC",
                      mergePeaks = TRUE, mapCorrect = TRUE, dataType = NULL){
  if(!is.null(dataType)){
    if(tolower(dataType) == "mnase" | tolower(dataType) == "sonicated"){
		#message("Using the default paraEM for MNase/sonicated data")
		minK=0;maxK=0;tol=1e-4;B=100;mSelect="AIC3";mergePeaks=TRUE;mapCorrect=TRUE;
    }
    else{
      stop("Invalid dataType")
    }
  }

  if(!is.finite(tol) & tol <= 0 & tol > 1){
    stop("'tol' must be a positive number between 0 and 1")
  }
  if(!is.finite(maxK) & maxK <= 0 & round(maxK) != maxK){
    stop("'maxK' must be a positive integer")
  }
  if(!is.finite(B) & B<1 & round(B) != B){
    stop("'B' must be a positive integer")
  }
  if(!is.character(mSelect) & ((mSelect!="BIC") | (mSelect!="AIC3"))){
    stop("'lambda' must be a positive number")
  }
  if(!is.logical(mergePeaks)){
    stop("'mergePeaks' must be a logical value")
  }
  if(!is.logical(mapCorrect)){
    stop("'mapCorrect' must be a logical value")
  }
  return(list(minK = minK, maxK = maxK, tol = tol, B = B, mSelect = mSelect,
              mergePeaks = mergePeaks, mapCorrect = mapCorrect))
}


#' Function that returns a list of parameters that can be used as an argument
#' of PICS.
#' 
#' This function takes from 0 to 6 parameters as argument, check if they are
#' valid and returns a list to be used in a call to PICS.
#' 
#' 
#' @param xi An \code{integer}. The average DNA fragment size.
#' @param rho An \code{integer}. A variance parameter for the average DNA
#' fragment size distribution.
#' @param alpha An \code{integer}. First hyperparameter of the inverse Gamma
#' distribution for sigma^2 in the PICS model
#' @param beta An \code{integer}. Second hyperparameter of the inverse Gamma
#' distribution for sigma^2 in the PICS model
#' @param lambda An \code{integer}. The precision of the prior for mu used for
#' histone data.
#' @param dMu An \code{integer}. Our best guess for the distance between two
#' neighboring nucleosomes.
#' @param dataType A character string. If a valid dataType is specified, use
#' our suggested parameters. "MNase" or "sonicated"
#' @param PExi A \code{numeric}. With paired end data, `xi' can be calculated
#' directly from the reads. If PExi is set, it will overwrite the xi determined
#' by the dataType.
#' @return Returns a list of 6 parameters to be used in PICS.
#' @author Renan Sauteraud
#' @seealso \code{PICS}
#' @keywords functions
#' @examples
#' 
#' # set prior for PICS data
#' paraPrior<-setParaPrior()
#' # set prior for sonicated data using our selected default parameters
#' paraPrior<-setParaPrior(dataType="sonicated")
#' 	
#' 
#' @export setParaPrior
setParaPrior <- function(xi = 200, rho = 1, alpha = 20, beta = 40000, lambda = 0,
                         dMu = 0, dataType = NULL, PExi = 0){
  if(!is.null(dataType)){
	  if(tolower(dataType) == "mnase"){
		  message("Using the default paraPrior for MNase data, for sonicated data use the argument dataType='sonicated'")
		  xi=150;rho=3;alpha=20;beta=20000;lambda=-0.000064;dMu=200; #Best default paramaeters
	  }
	  else if(tolower(dataType)=="sonicated"){
		  message("Using the default paraPrior for sonicated data")
		  xi=150;rho=1.2;alpha=10;beta=20000;lambda=-0.000064;dMu=200; #Best default parameters
	  }
	  else{
		  stop("Invalid dataType")
	  }
  }
  if(!is.finite(xi)){
    stop("'xi' must be a numeric value")
  }
  if(!is.finite(rho) & rho<=0){
    stop("'rho' must be a positive number")
  }
  if(!is.finite(alpha) & alpha<=0){
    stop("'alpha' must be a positive number")
  }
  if(!is.finite(beta) & beta<=0){
    stop("'beta' must be a positive number")
  }
  if(!is.finite(lambda) & lambda<=0){
    stop("'lambda' must be a positive number")
  }
  if(!is.finite(dMu) & dMu<=0){
    stop("'dMu' must be a positive number")
  }
  if(PExi>0){
    xi<-PExi
  }
  return(list(xi = xi, rho = rho, alpha = alpha, beta = beta,
              lambda = lambda, dMu = dMu))
}