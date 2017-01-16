######################
#clear console       #
######################
#rm(list=ls())

######################
#load  R libraries   #
######################
packages <- list(quadprog = "quadprog",
                 quantreg = "quantreg",
                 VGAM = "VGAM",
                 MCMCpack = "MCMCpack",
                 splines = "splines",
                 ald = "ald",
                 emulator = "emulator")
lapply(packages, require, character.only = T)
# List of functions:
## 1) checkCharacter
## 2) checkBinary
## 3) checkNumericValues
## 4) CheckScalar
## 6) CheckNonnegativeNumber
## 7) CheckPositiveNumber
## 8) CheckPositiveScalar
## 9) checkInteger
## 10) CheckUnitInverval
## 11) checkLocs
## 12) checkInteriorKnots
## 13) CheckInvertibileSymmetricMatrix

## 14) setEqual
## 14) calculateMeanan
## 15) calculateVariance
## 16) findIntegerMax
## 16) findInterval
## 17) matrixMultiply
## 18) calculateKronecker
## 19) quadraticForm
## 20) quadraticFormcalculateKronecker
## 21) cholMatrix
## 22) invertSymmetricMatrix
## 23) makeAbsoluteDistance
## 24) makeAutoregressivePrecision
## 25) makeSpatialDistance
## 26) makeSpatialPrecision

## 27) makeSplineKnots
## 28) mSpline
## 29) iSpline
## 30) cubicMSpline
## 31) cubicISpline
## 32) makeCubicSplineCoefficients
## 33) cubicMSpline2
## 34) cubicISpline2
## 35) makeStickBreakingWeights

## 36) logDensityNormal
## 37) logDensityT
## 38) logDLogistic
## 39) LogDAsymmetricLapace
## 40) logDensityWeibull
## 41) logDensityGamma
## 42) cdfNormal
## 43) cdfT
## 44) cdfLogistic
## 45) PAsymmetricLapace
## 46) cdfWeibull
## 47) pGamma
## 48) qNorm
## 49) qT
## 50) qLogistic
## 51) QAsymmetricLapace
## 52) qWeibull
## 53) qGamma
## 54) truncatedGaussianMoments

## 55) randomNormal
## 56) randomTruncatedNormal
## 57) randomInverseWishart
## 58) randomDirichlet
## 59) randomMultinomial
## 60) randomBeta
## 61) randomGamma
## 62) UpdatecalculateKroneckerInverseWishart

checkCharacter <- function(v, vName = NULL) {
  # Checks that an argument is character.
  #
  # Args:
  #   v: an object to be tested if binary.
  #   vName: the name of v passed along from another function.
  # Returns nothing.
  #
  # Error handling:
  vName <- ifelse(is.null(vName), deparse(substitute(v)), vName)
  if(!is.character(v)) {
    stop(message =
           paste0(vName, " must be a character.\n"))
  }
}

checkBinary <- function(v, vName = NULL, noMissing = T, lengthOne = T) {
  # Checks that an argument is binary 0/1.
  #
  # Args:
  #   v: an object to be tested if binary.
  #   vName: the name of v passed along from another function.
  #   noMissing: an indicator that no missing values are allowed.
  #   lengthOne: an indicator that the argument is of length 1.
  # Returns nothing.
  #
  # Error handling:
  if(length(noMissing) != 1) {
    stop(message = "noMissing must be of length 1.\n")
  }
  check1 <-  noMissing %in% c(T, F)
  if(!check1) {
    stop(message = "noMissing must be either TRUE or FALSE.\n")
  }
  if(length(lengthOne) != 1) {
    stop(message = "lengthOne must be of length 1.\n")
  }
  check2 <-  lengthOne %in% c(T, F)
  if(!check2) {
    stop(message = "lengthOne must be either TRUE or FALSE.\n")
  }
  vName <- ifelse(is.null(vName), deparse(substitute(v)), vName)
  checkCharacter(vName)
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(vName, " cannot contain only NA.\n"))
  }
  if(noMissing) {
    if(any(is.na(v))) {
      stop(message =
             paste0(vName, " must not have missing values.\n"))
    }
  }
  if(lengthOne) {
    if(length(v) != 1) {
      stop(message = paste0(vName, " must be of length 1.\n"))
    }  
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
  zeroTest <- v == 0
  oneTest <- v == 1
  if(any(!zeroTest & !oneTest)) {
    stop(message = paste0(vName, " must be composed of 0 or 1.\n"))
  }
}

checkInteger <- function(v, vName = NULL, noMissing = T, lengthOne = T,
                         nonnegative = F, positive = F) {
  # Checks that a scalar argument is an integer.
  #  It does not have to be of integer type.
  # Args:
  #   v: an object to be tested if numeric.
  #   vName: the name of v passed along from another function.
  #   noMissing: an indicator that no missing values are allowed.
  #   lengthOne: an indicator that the argument is of length 1.
  #   nonnegative: an indicator that values should be nonnegative.
  #   positive: an indicator that values should be positive.
  # Returns nothing.
  #
  # Error handling:
  vName <- ifelse(is.null(vName), deparse(substitute(v)), vName)
  checkCharacter(vName)
  checkBinary(noMissing)
  checkBinary(lengthOne)
  checkBinary(nonnegative)
  checkBinary(positive)
 
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(vName, " cannot contain only NA.\n"))
  }
  if(noMissing) {
    if(any(is.na(v))) {
      stop(message =
             paste0(vName, " must not have missing values.\n"))
    }
  }
  if(lengthOne) {
    if(length(v) != 1) {
      stop(message = paste0(vName, " must be of length 1.\n"))
    }  
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
  if(any(v %% 1 != 0)) {
    stop(message =
           paste0(vName, " must be an integer.\n"))
  }
  if(nonnegative) {
    if(any(v < 0)) {
      stop(message =
             paste0(vName, " must have nonnegative elements.\n"))
    }
  }
  if(positive) {
    if(any(v <= 0)) {
      stop(message =
             paste0(vName, " must have positive elements.\n"))
    }
  }
}

checkNumeric <- function(v, vName = NULL, noMissing = T, lengthOne = T,
                         nonnegative = F, positive = F) {
  # Checks that an argument is numeric.
  #
  # Args:
  #   v: an object to be tested if numeric.
  #   vName: the name of v passed along from another function.
  #   noMissing: an indicator that no missing values are allowed.
  #   lengthOne: an indicator that the argument is of length 1.
  #   nonnegative: an indicator that values should be nonnegative.
  #   positive: an indicator that values should be positive.
  # Returns nothing.
  #
  # Error handling:
  vName <- ifelse(is.null(vName), deparse(substitute(v)), vName)
  checkCharacter(vName)
  checkBinary(noMissing)
  checkBinary(lengthOne)
  checkBinary(nonnegative)
  checkBinary(positive)
 
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(vName, " cannot contain only NA.\n"))
  }
  if(noMissing) {
    if(any(is.na(v))) {
      stop(message =
             paste0(vName, " must not have missing values.\n"))
    }
  }
  if(lengthOne) {
    if(length(v) != 1) {
      stop(message = paste0(vName, " must be of length 1.\n"))
    }  
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
 
  if(!is.numeric(v)) {
    stop(message =
           paste0(vName, " must be numeric.\n"))
  }
  if(nonnegative) {
    if(any(v < 0)) {
      stop(message =
             paste0(vName, " must have nonnegative elements.\n"))
    }
  }
  if(positive) {
    if(any(v <= 0)) {
      stop(message =
             paste0(vName, " must have positive elements.\n"))
    }
  }
}

checkUnitInterval <- function(v, vName = NULL,
                              noMissing = T, lengthOne = T, equal = F) {
  # function that checks that a number is in (0,1) or [0,1].
  #
  # Args:
  #   v: an object to be tested if numeric.
  #   vName: the name of v passed along from another function.
  #   noMissing: an indicator that no missing values are allowed.
  #   lengthOne: an indicator that the argument is of length 1.
  #   equal: an indicator that the valid interval is [0, 1], not (0, 1).
  # Returns nothing.
  #
  # Error handling:
  vName <- ifelse(is.null(vName), deparse(substitute(v)), vName)
  checkCharacter(vName)
  checkBinary(noMissing)
  checkBinary(lengthOne)
  checkBinary(equal)
 
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(vName, " cannot contain only NA.\n"))
  }
  if(noMissing) {
    if(any(is.na(v))) {
      stop(message =
             paste0(vName, " must not have missing values.\n"))
    }
  }
  if(lengthOne) {
    if(length(v) != 1) {
      stop(message = paste0(vName, " must be of length 1.\n"))
    }  
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
  if(!equal) {
    if(any(v <= 0)) {
      stop(message = paste0(vName, " cannot be less than or equal to 0.\n"))
    }
    if(any(v >= 1)) {
      stop(message = paste0(vName, " cannot be greater than or equal to 1.\n"))
    }  
  } else {
    if(any(v < 0)) {
      stop(message = paste0(vName, " cannot be less than 0.\n"))
    }
    if(any(v > 1)) {
      stop(message = paste0(vName, " cannot be greater than 1.\n"))
    }  
  }
}

checkLocs <- function(locs) {
  # Check that locs are valid.  
  # Args:
  #   locs: matrix of locations.
  # Returns nothing.
  #
  if(!is.matrix(locs)) {
    stop(message = "locs must be a matrix.\n")
  }
  if(ncol(locs) != 2) {
    stop(message = "locs must have 2 columns.\n")
  }
  checkNumeric(locs, lengthOne = F)
  if(length(locs) != length(unique(locs))) {
    stop(message = "locs must be unique.\n")
  }
}

checkInteriorKnots <- function(interiorKnots) {
  # check that interior knot sequence is of correct form.
  # Args:
  #   interiorKnots: vector of knots that are in (0, 1).
  # Returns nothing
  # Error handling
  checkUnitInterval(interiorKnots, lengthOne = F)
  nKnots <- length(interiorKnots)
  if(sum(order(interiorKnots) == (1:nKnots)) < nKnots) {
    stop(message = "interiorKnots must be ordered sequentially.\n")
  }
}

checkCubicSplineKnots <- function(splineKnots) {
  # check that cubic spline knots is of correct form.
  # Args:
  #   splineKnots: vector of knots.
  # Returns nothing
  # Error handling
  checkUnitInterval(splineKnots, lengthOne = F, equal = T)
  nKnots <- length(splineKnots)
  if(nKnots < 9) {
    stop(message = "splineKnots must be a sequence of length at least 9.\n")
  }
  if(any(splineKnots[1:3] != 0)) {
    stop(message = "splineKnots must have first 3 elements equal to 0.\n")
  }
  if(any(splineKnots[(nKnots - 4):nKnots] != 1)) {
    stop(message = "splineKnots must have last 5 elements equal to 1.\n")
  }
  interiorKnots <- splineKnots[4:(nKnots - 5)]
  checkUnitInterval(interiorKnots, lengthOne = F)
  if(any(order(interiorKnots) != 1:length(interiorKnots))) {
    stop(message = "splineKnots must be sorted in increasing order.\n")
  }
}

checkInvertibleSymmetricMatrix <- function(v, vName = NULL) {
  # Checks that a matrix is invertible and symmetric.
  #
  # Args:
  #   v: matrix.
  # Returns nothing.
  # Error handling:
  vName <- ifelse(is.null(vName), deparse(substitute(v)), vName)
  checkCharacter(vName)
  if(!is.matrix(v)) {
    stop(message =
           paste0(vName, " must be a matrix.\n"))
  }
  if(!isSymmetric(v)) {
    stop(message =
           paste0(vName, " must be symmetric.\n"))
  }
  flag <- 0
  vInv <- 0
  suppressWarnings(
    tryCatch(
      vInv <- chol2inv(
        chol(v)
      ),
      error = function(e) {return(vInv)}
    )
  )
  if (length(vInv) == 1) {
    if (vInv == 0) {
      flag <- 1
    }
  }
  if(flag == 1) {
    stop(message =
           paste0(vName, " is symmetric but not invertible.\n"))
  }
}

setEqual <- function(x, cIndicator = F) {
# creates a copy of numeric array x.
  # Args:
  #    x: vector of reals.
  # Returns a copy of x.
  # Error handling:
  checkNumeric(x, lengthOne = F)
  checkBinary(cIndicator)
  if(!cIndicator) {
    return(x)
  } else {
    lengthX <- length(x)
    out <- .C("setEqual",
              lengthX = as.integer(lengthX),
              x = as.double(x),
              copyX  = as.double(rep(0, lengthX))
    )
    x <- out$x
  }
  return(x)
}

prettyPlot <- function(x, y) {
  # prettyPlots vector y vs vector x.
  # Args:
  #    x: vector of reals.
  #    y: vector of reals.
  # Returns nothing.
  # Error handling:
  checkNumeric(x, lengthOne = F)
  checkNumeric(y, lengthOne = F)
  if(length(x) != length(y)) {
    stop(message = "length of x must equal length of y.\n")
  }
  # mar: numeric vector of length 4, which sets the margin sizes
  ## in the following order: bottom, left, top, and right.
  ## Default is c(5.1, 4.1, 4.1, 2.1).
  # mgp: A numeric vector of length 3, which sets the axis label locations
  ## relative to the edge of the inner plot window.
  ## The first value represents the location the labels (i.e. xlab and ylab in plot),
  ## the second the tick-mark labels, and third the tick marks.
  ## The default is c(3, 1, 0).
  # oma: outer margin area.  Provides space around inner plot area.  Default is rep(0, 4).
  par(mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0), tck = -0.01)
  plot(x, y)
}

calculateMean <- function(x, cIndicator = F) {
  # Returns the mean of vector x.
  # Args:
  #    x: vector of reals.
  #    cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   meanX: the mean of x.
  # Error handling:
  checkNumeric(x, lengthOne = F)
  checkBinary(cIndicator)
  if(!cIndicator) {
    meanX <- mean(x)
  } else {
    out <- .C("calculateMean",
              lengthX = as.integer(length(x)),
              x = as.double(x),
              meanX  = as.double(0)
    )
    meanX <- out$meanX
  }
  return(meanX)
}

calculateVariance <- function(x, cIndicator = F) {
  # Returns the variance of vector x.
  # Args:
  #    x: vector of reals.
  #    cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   meanX: the mean of x.
  # Error handling:
  checkNumeric(x, lengthOne = F)
  checkBinary(cIndicator)
  lengthX <- length(x)
  if(lengthX == 1) {
    stop(message = "x must be of length greater than 1.\n")
  }
  if(!cIndicator) {
    varX <- var(x)
  } else {
    out <- .C("calculateVariance",
              lengthX = as.integer(length(x)),
              x = as.double(x),
              varX  = as.double(0)
    )
    varX <- out$varX
  }
  return(varX)
}

findIntegerMax <- function(x, cIndicator = F) {
  # finds the maximum value in a vector of integers.
  # Args:
  #   x: the vector of integers to search through.
  # Returns:
  #   maxX: the maximum value of x.
  # Error handling:
  checkInteger(x)
  checkBinary(cIndicator)
  if(!cIndicator) {
    maxX <- max(x)
  } else {
    lengthX <- length(x)
    out <- .C("findIntegerMax",
              lengthX = as.integer(lengthX),
              x = as.integer(x),
              maxX = as.integer(0)
    )
    maxX <- out$maxX
  }
  return(maxX)
}

findInterval <- function(tau, x, cIndicator = F) {
  # Finds which interval of x scalar tau is in
  #
  # Args:
  #   tau: real scalar
  #   x: sorted vector of reals
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   bin: the bin (indexed from 0 to length(x) - 1) such that x[bin] <= tau < x[bin + 1].
  # For example: if x = 1:10, and tau = 1.5, then bin = 0;
  #              if x = 1:10, and tau = 9.5, then bin = 8;
  #              if x = 1:10, and tau = 100, then bin = 9;
  #  In R add 1 to the bins above.
  # Error handling:
  checkNumeric(tau)
  checkNumeric(x, lengthOne = F)
  lengthX <- length(x)
  if(lengthX != length(unique(x))) {
    stop(message = "values of x must be unique.\n")
  }
  if(any(order(x) != 1:lengthX)) {
    stop(message = "x must be in ascending order.\n")
  }
 
  if(tau < x[1]) {
    stop(message = "tau must be greater than or equal to lowest value of x.\n")
  }
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    bin <- findInterval(tau, x)
  } else {
    out <- .C("findInterval",
              lengthX = as.integer(lengthX),
              x = as.double(x),
              tau  = as.double(tau),
              bin = as.integer(1)
    )
    bin <- out$bin
  }
  return(bin)
}

matrixMultiply <- function(a, b, transA = 0, transB = 0, cIndicator = F) {
  # Calculates the product of matrices a and b.
  #
  # Args:
  #   a: matrix of reals.
  #   b: matrix of reals.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   c: product of a and b.
  #
  # Error handling:
  checkNumeric(a, lengthOne = F)
  checkNumeric(b, lengthOne = F)
  if(!is.matrix(a)) {
    stop(message = "a must be a matrix.\n")
  }
  if(!is.matrix(b)) {
    stop(message = "b must be a matrix.\n")
  }
  nRowA <- nrow(a)
  nColA <- ncol(a)
  nRowB <- nrow(b)
  nColB <- ncol(b)
 
  checkBinary(transA)
  checkBinary(transB)
 
  if(transA == 0 && transB == 0) {
    if(nColA != nRowB) {
      stop(message = "number of columns of a must equal number of rows of b.\n")
    }
  }
  if(transA == 1 && transB == 0) {
    if(nRowA != nRowB) {
      stop(message = "number of rows of a must equal number of rows of b.\n")
    }
  }
  if(transA == 0 && transB == 1) {
    if(nColA != nColB) {
      stop(message = "number of columns of a must equal number of cols of b.\n")
    }
  }
  if(transA == 1 && transB == 1) {
    if(nRowA != nColB) {
      stop(message = "number of rows of a must equal number of cols of b.\n")
    }
  }
 
  checkBinary(cIndicator)
  if(!cIndicator) {
    if(transA) {
      a <- t(a)
    }
    if(transB) {
      b <- t(b)
    }
    c <- a %*% b
  } else {
    nRowA <- ifelse(transA, nColA, nRowA)
    nColC <- ifelse(transB, nRowB, nColB)
    midC <- ifelse(transA, nRowA, nColA)
    dimC <- c(nRowA, nColC)
    out <- .C("matrixMultiply",
              transA = as.integer(transA),
              transB = as.integer(transB),
              nRowA  = as.integer(nRowA),
              nColA  = as.integer(nColA),
              nRowB  = as.integer(nRowB),
              nColB  = as.integer(nColB),
              a      = as.double(a),
              b      = as.double(b),
              c      = as.double(rep(0, nRowA * nColC))
    )
    c <- array(out$c, dim = dimC)
  }  
  return(c)
}

calculateKronecker <- function(a, b, cIndicator = F) {
  # Calculates the kronecker product of a and b.
  #
  # Args:
  #   a: matrix of reals.
  #   b: matrix of reals.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   c: kronecker product of a and b.
  #
  # Error handling:
  checkNumeric(a, lengthOne = F)
  checkNumeric(b, lengthOne = F)
  if(!is.matrix(a)) {
    stop(message = "a must be a matrix.\n")
  }
  if(!is.matrix(b)) {
    stop(message = "b must be a matrix.\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    c <- kronecker(a, b)
  } else {
    nRowA <- nrow(a)
    nColA <- ncol(a)
    nRowB <- nrow(b)
    nColB <- ncol(b)
   
    out <- .C("calculateKronecker",
              nRowA  = as.integer(nRowA),
              nColA  = as.integer(nColA),
              nRowB  = as.integer(nRowB),
              nColB  = as.integer(nColB),
              a       = as.double(a),
              b       = as.double(b),
              c       = as.double(rep(0, nRowA * nRowB * nColA * nColB))
    )  
    c <- array(out$c, dim = c(nRowA * nRowB, nColA * nColB))
  }
  return(c)
}

quadraticForm <- function(x, a, cIndicator = F) {
  # Finds qf = x'a x, where a is symmetric.
  #
  # Args:
  #   x: vector.
  #   a: symmetric matrix that is length(x) x length(x).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   qf: the quadratic form.
  # Error handling:
  checkNumeric(x, lengthOne = F)
  lengthX <- length(x)
  checkNumeric(a, lengthOne = F)
  if(!is.matrix(a)) {
    stop(message = "a must be a matrix.\n")
  }
  if(!isSymmetric(a)) {
    stop(message = "a must be symmetric.\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    qf <- t(x) %*% a %*% x
  } else {
    out <- .C("quadraticForm",
              lengthX = as.integer(lengthX),
              x = as.double(x),
              a  = as.double(a),
              qf = as.double(0)
    )
    qf <- out$qf
  }
  return(qf)
}

quadraticFormcalculateKronecker <- function(dMat, eMat, x, cIndicator = F) {
  # This function returns t(x) %*% kronecker(dMat, eMat) %*% x
  #   without calculating the kronecker product
  #
  # Args:
  #   dMat : square matrix.
  #   eMat : square matrix.
  #   x: vector of length nrow(dMat) * nrow(eMat).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns the quadratic form.
  # Error handling:
  checkNumeric(dMat, lengthOne = F)
  if(!is.matrix(dMat)) {
    stop(message = "dMat must be a matrix.\n")
  }
  checkNumeric(eMat, lengthOne = F)
  if(!is.matrix(eMat)) {
    stop(message = "eMat must be a matrix.\n")
  }
  nRowD <- nrow(dMat)
  if(nRowD != ncol(dMat)) {
    stop(message = "dMat must be a square matrix.\n")
  }
  nRowE <- nrow(eMat)
  if(nRowE != ncol(eMat)) {
    stop(message = "eMat must be a square matrix.\n")
  }
  checkNumeric(x, lengthOne = F)
  if(length(x) != (nRowD * nRowE)) {
    stop(message = "x must be of length nrow(d) * nrow(e).\n")
  }
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    newX <- matrix(x, nRowE, nRowD)
    matrixProduct <- eMat %*% newX %*% dMat
    qf <- sum(newX * matrixProduct)
   
    #     qf <- 0
    #     for(d in 1:nRowD) {
    #       qf <- qf + newX[, d] %*% matrixProduct[, d]
    #     }
  } else {
    out <- .C("quadraticFormcalculateKronecker",
              nRowD  = as.integer(nRowD),
              nRowE  = as.integer(nRowE),
              d       = as.double(dMat),
              e       = as.double(eMat),
              x       = as.double(x),
              qf      = as.double(0)
             
    )
    qf <- out$qf
  }
  return(qf)
}

cholMatrix <- function(v, cIndicator = F) {
  # Calculates the cholesky decompostion of v and the log determinant of v
  #
  # Args:
  #   v: a positive definite symmetric matrix
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   vChol: the lower triangular Cholesky decomposition of v.
  #   logDet: the log determinant of v.
  # Some error handling done in C.
  checkNumeric(v, lengthOne = F)
  if(!is.matrix(v)) {
    stop(message = "v must be a matrix.\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    vChol <- t(chol(v))
    logDet <- 2 * sum(log(diag(vChol)))
  }
  else {
    nRowV <- nrow(v)
    out <- .C("cholMatrix",
              nRowV   = as.integer(nRowV),
              v       = as.double(v),          
              vChol   = as.double(rep(0, nRowV * nRowV)),
              logDet  = as.double(0)
    )
    vChol <- array(out$vChol, dim = c(nRowV, nRowV))
    logDet <- out$logDet
  }
  return(
    list(vChol = vChol, logDet = logDet)
  )
}

invertSymmetricMatrix <- function(v, cIndicator = F) {
  # Calculates the inverse of v and the log determinant of the inverse.
  #
  # Args:
  #   v: an invertible symmetric matrix.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   vInv: the inverse of v.
  #   logDetVInv: the log determinant of the inverse.
  # Some error handling done in C.
  checkNumeric(v, lengthOne = F)
  if(!is.matrix(v)) {
    stop(message = "v must be a matrix.\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    cholV <- chol(v) #upper triangular form
    logDet <- -2 * sum(log(diag(cholV)))
    vInv <- chol2inv(cholV)
  } else {
    nRowV <- nrow(v)
    out <- .C("invertSymmetricMatrix",
              nRowV   = as.integer(nRowV),
              v       = as.double(v),          
              vInv    = as.double(rep(0, nRowV * nRowV)),
              logDet  = as.double(0)
    )
    vInv <- array(out$vInv, dim = c(nRowV, nRowV))
    logDet <- out$logDet  
  }
  list(vInv = vInv, logDet = logDet)
}

makeAbsoluteDistance <- function(timepoints, cIndicator = F) {
  # Constructs absolute distance matrix using unique vector of timepoints.
  #
  # Args:
  #   timepoints: vector of unique timepoints.
  # Returns:
  #   distanceMatrix: autoregressive distance matrix.
  #   Error handling:
  if(FALSE %in% (timepoints == sort(timepoints))) {
    stop(message = "timepoints must be ordered sequentially.\n")
  }
  if(length(timepoints) != length(unique(timepoints))) {
    stop(message = "timepoints must contain unique values.\n")
  }
  if(!cIndicator) {
    distanceMatrix <- abs(outer(timepoints, timepoints, "-"))
  } else {
    nTimepoints <- length(timepoints)
    out <- .C("makeAbsoluteDistance",
              nTimepoints = as.integer(nTimepoints),
              timepoints   = as.double(timepoints),
              distanceMatrix = as.double(rep(0, nTimepoints * nTimepoints))
    )
    distanceMatrix <- matrix(out$distanceMatrix, nTimepoints, nTimepoints)
  }
  return(distanceMatrix)
}

makeAutoregressivePrecision <- function(rho, distanceMatrix, cIndicator = F){
  # Constructs precision matrix of AR-1 process.
  #
  # Args:
  #   rho: scalar in (0, 1).  
  #   distanceMatrix: matrix indicating the distance between timepoints.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   precAR: autoregressive precision matrix.
  #   logDet: log determinant of the precision matrix.
  # Error handling:
  if(length(rho) > 1){
    stop(message = "rho must be of length 1.\n")
  }
  checkUnitInterval(rho)
  if(!is.matrix(distanceMatrix)) {
    stop(message = "distanceMatrix must be a matrix.\n")
  }
  checkNumeric(distanceMatrix, lengthOne = F)
  if(!isSymmetric(distanceMatrix)) {
    stop(message = "distanceMatrix must be symmetric.\n")
  }
  if(any(diag(distanceMatrix) != 0)) {
    stop(message = "distanceMatrix must have a diagonal of 0.\n")
  }
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    ar <- rho ^ distanceMatrix
    out <- invertSymmetricMatrix(ar)
    precAR <- out$vInv
    logDet <- out$logDet
  } else {
    nTimepoints <- nrow(distanceMatrix)
    out <- .C("makeAutoregressivePrecision",
              nTimepoints  = as.integer(nTimepoints),
              rho           = as.double(rho),
              dist          = as.double(distanceMatrix),
              precAR       = as.double(rep(0, nTimepoints * nTimepoints)),
              logDet       = as.double(0)
    )
    precAR <- array(out$precAR, dim = c(nTimepoints, nTimepoints))
    logDet <- out$logDet
  }
  return(list(precAR = precAR, logDet = logDet))
}

makeSpatialDistance <- function(locs) {
  # Constructs Euclidean distance matrix using unique vector of locations.
  #
  # Args:
  #   locs: matrix with 2 columns of unique locations.
  # Returns:
  #   spatialDist: Euclidean distance matrix.
  # Error handling:
  checkLocs(locs)
  nLocs <- nrow(locs)
  distLon <- outer(locs[,1], locs[,1], "-")
  distLat <- outer(locs[,2], locs[,2], "-")
  spatialDist <- sqrt(distLon ^ 2 + distLat ^ 2)
  return(spatialDist)
}

makeSpatialPrecision <-function(rho, spatialDist, cIndicator = F){
  # Constructs exponential spatial precision matrix.
  #
  # Args:
  #   rho: range parameter (positive scalar).
  #   spatialDist: euclidean distance matrix.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   precSpatial: spatial precision matrix.
  #   logDet: log determinant of the precision matrix.
  #
  # Error handling:
  if(length(rho) > 1){
    stop(message = "rho must be of length 1.\n")
  }
  checkNumeric(rho, lengthOne = T, positive = T)
  if(!is.matrix(spatialDist)) {
    stop(message = "spatialDist must be a matrix.\n")
  }
  checkNumeric(spatialDist, lengthOne = F)
  if(!isSymmetric(spatialDist)) {
    stop(message = "spatialDist must be symmetric.\n")
  }
  if(any(diag(spatialDist) != 0)) {
    stop(message = "spatialDist must have a diagonal of 0.\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    sp <- exp(-spatialDist / rho)
    out <- invertSymmetricMatrix(sp)
    precSpatial <- out$vInv
    logDet <- out$logDet
  } else {
    nLocs <- nrow(spatialDist)
    out <- .C("makeSpatialPrecision",
              nLocs       = as.integer(nLocs),
              rho         = as.double(rho),
              dist        = as.double(spatialDist),
              precSpatial = as.double(rep(0, nLocs * nLocs)),
              logDet      = as.double(0)
    )
    precSpatial <- array(out$precSpatial, dim = c(nLocs, nLocs))
    logDet <- out$logDet
  }
  return(list(precSpatial = precSpatial, logDet = logDet))
}

makeSplineKnots <- function(interiorKnots, splineDegreesOfFreedom = 3, cIndicator = F) {
  #constructs knot sequence for spline basis.
  # Args:
  #   interiorKnots = vector of knots that are in (0, 1).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   splineKnots: vector of knots for the spline basis.
  # Error handling:
  checkInteriorKnots(interiorKnots)
  checkInteger(splineDegreesOfFreedom, positive = T)
  checkBinary(cIndicator)
  if(!cIndicator) {
    splineKnots <- c(rep(0, splineDegreesOfFreedom), interiorKnots, rep(1, splineDegreesOfFreedom + 2))  
  } else {
    interiorKnotsLength <- length(interiorKnots)
    out <- .C("makeSplineKnots",
              interiorKnotsLength     = as.integer(interiorKnotsLength),
              splineDegreesOfFreedom  = as.integer(splineDegreesOfFreedom),
              interiorKnots           = as.double(interiorKnots),
              splineKnots             = as.double(rep(0, 2 * splineDegreesOfFreedom + interiorKnotsLength + 2))
    )
    splineKnots <- out$splineKnots
  }  
  return(splineKnots)
}

mSpline <- function(tau, splineDegreesOfFreedom, splineKnots, m, cIndicator = F) {
  # Given knots splineKnots, evaluates the mth m-spline of
  #   splineDegreesOfFreedom degrees of freedom at tau.
  # Args:
  #   tau: scalar in (0, 1).
  #   splineDegreesOfFreedom: degrees of freedom for splines.
  #   splineKnots: knot sequence.
  #   m: index for the spline.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interiorKnots and splineDegreesOfFreedom degrees of freedom,
  #       the mth m-spline evaluated at tau.
  # Error handling
  checkUnitInterval(tau)
  checkInteger(splineDegreesOfFreedom, positive = T)
  checkNumeric(splineKnots, lengthOne = F, nonnegative = T)
  if(min(splineKnots) != 0) {
    stop(message = "splineKnots first element must be 0.\n")
  }
  if(max(splineKnots) != 1) {
    stop(message = "splineKnots greatest element must be 1.\n")
  }
  if(any(order(splineKnots) != (1:length(splineKnots)))) {
    stop(message = "splineKnots must be in increasing order.\n")
  }
 
  checkInteger(m, positive = T)
  if(m  > length(splineKnots) - splineDegreesOfFreedom + 1) {
    stop(message = "m cannot be greater than length(splineKnots) - splineDegreesOfFreedom + 1.\n")
  }
 
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    if(tau < splineKnots[m] || tau >= splineKnots[m + splineDegreesOfFreedom]) {
      v <- 0
    } else if (splineDegreesOfFreedom == 1) {
      v <- 1 / (splineKnots[m + 1] - splineKnots[m])
    } else {
      d1 <- tau - splineKnots[m]
      d2 <- splineKnots[m + splineDegreesOfFreedom] - tau
      v <- d1 * mSpline(tau, splineDegreesOfFreedom - 1, splineKnots, m) +
        d2 * mSpline(tau, splineDegreesOfFreedom - 1, splineKnots, m + 1)
      v <- v * splineDegreesOfFreedom
      v <- v / ((splineDegreesOfFreedom - 1) * (splineKnots[m + splineDegreesOfFreedom] - splineKnots[m]))
    }
  } else {
    out <- .C("mSpline",
              tau                     = as.double(tau),          
              splineDegreesOfFreedom  = as.integer(splineDegreesOfFreedom),          
              splineKnots             = as.double(splineKnots),          
              m                       = as.integer(m - 1),
              v                       = as.double(0)          
    )  
    v <- out$v
  }
  return(v)
}

iSpline <- function(tau, splineDegreesOfFreedom, splineKnots, m, cIndicator = F) {
  # Given knots splineKnots, evaluates the mth i-spline of
  #   splineDegreesOfFreedom degrees of freedom at tau.
  # Args:
  #   tau: scalar in (0, 1).
  #   splineDegreesOfFreedom: degrees of freedom for splines.
  #   splineKnots: knot sequence.
  #   m: index for the spline.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interiorKnots and splineDegreesOfFreedom degrees of freedom,
  #       the mth i-spline evaluated at tau.
  # Error handling
  checkUnitInterval(tau)
  checkInteger(splineDegreesOfFreedom, positive = T)
 
  checkNumeric(splineKnots, lengthOne = F, nonnegative = T)
  if(min(splineKnots) != 0) {
    stop(message = "splineKnots first element must be 0.\n")
  }
  if(max(splineKnots) != 1) {
    stop(message = "splineKnots greatest element must be 1.\n")
  }
  if(any(order(splineKnots) != (1:length(splineKnots)))) {
    stop(message = "splineKnots must be in increasing order.\n")
  }
 
  checkInteger(m, positive = T)
  if(m  > length(splineKnots) - splineDegreesOfFreedom) {
    stop(message = "m cannot be greater than length(splineKnots) - 3.\n")
  }
 
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    bin <- findInterval(tau, splineKnots)
    if (bin < m){
      v<- 0
    } else if ((bin - splineDegreesOfFreedom + 1) > m ){
      v <- 1
    } else {
      v <- 0
      for (l in m : (bin + 1)){
        v <- v + mSpline(tau, splineDegreesOfFreedom + 1, splineKnots, l) *
          (splineKnots[l + splineDegreesOfFreedom + 1] - splineKnots[l]) /
          (splineDegreesOfFreedom + 1)
      }
    }
  } else {
    splineKnotsLength <- length(splineKnots)
    out <- .C("iSpline",
              tau                     = as.double(tau),          
              splineDegreesOfFreedom  = as.integer(splineDegreesOfFreedom),          
              splineKnotsLength       = as.integer(splineKnotsLength),          
              splineKnots             = as.double(splineKnots),          
              m                       = as.integer(m - 1),
              v                       = as.double(0)          
    )  
    v <- out$v
  }
  return(v)
}

cubicMSpline <- function(tau, interiorKnots, m, cIndicator = F) {
  # Given knots splineKnots, evaluates the mth m-spline of
  #   3 degrees of freedom at tau.
  # Args:
  #   tau: scalar in (0, 1).
  #   splineDegreesOfFreedom: degrees of freedom for splines (should be 3).
  #   interiorKnots: interior knot sequence.
  #   m: index for the spline.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interiorKnots and splineDegreesOfFreedom degrees of freedom,
  #       the mth m-spline evaluated at tau
  # Error handling
  checkUnitInterval(tau)
  checkInteriorKnots(interiorKnots)
  checkInteger(m, positive = T)
  splineKnots <- makeSplineKnots(interiorKnots, splineDegreesOfFreedom = 3, cIndicator = F)
 
  if(m  > length(splineKnots) - 3) {
    stop(message = "m cannot be greater than length(splineKnots) - 3.\n")
  }
  checkBinary(cIndicator)
  splineDegreesOfFreedom <- 3
  if(!cIndicator) {
    bin <- findInterval(tau, splineKnots)
    if (bin < m) {
      v <- 0
    } else if (
      (bin - splineDegreesOfFreedom + 1) > m) {
      v <- 0
    } else if (bin == m) {
      v <- 3 * (tau - splineKnots[m]) ^ 2 /
        ((splineKnots[m + 1] - splineKnots[m]) *
           (splineKnots[m + 2] - splineKnots[m]) *
           (splineKnots[m + 3] - splineKnots[m])
        )
    } else if (bin == (m + 1)) {
      i1 <-  1 / (splineKnots[m + 2] - splineKnots[m + 1])
      i2 <-  (2 * (splineKnots[m + 3] - splineKnots[m + 1])) /
        ((splineKnots[m + 3] - splineKnots[m + 1]) *
           (splineKnots[m + 2] - splineKnots[m + 1])
        )
      i3 <-  -3 *(
        (splineKnots[m + 3] - tau) ^ 2/
          ((splineKnots[m + 3] - splineKnots[m + 1]) *
             (splineKnots[m + 3] - splineKnots[m]) *
             (splineKnots[m + 2] - splineKnots[m + 1])) +
          (tau - splineKnots[m]) ^ 2/
          (
            (splineKnots[m + 3] - splineKnots[m]) *
              (splineKnots[m + 2] - splineKnots[m]) *
              (splineKnots[m + 2] - splineKnots[m + 1])
          )    
      )
      v <- i1 + i2 + i3
    } else{
      v <- 3 * (splineKnots[m + 3] - tau) ^ 2 /
        ((splineKnots[m + 3] - splineKnots[m + 2]) *
           (splineKnots[m + 3] - splineKnots[m + 1]) *
           (splineKnots[m + 3] - splineKnots[m])
        )
    }
  } else {
    out <- .C("cubicMSpline",
              tau                 = as.double(tau),          
              splineKnotsLength   = as.integer(length(splineKnots)),          
              splineKnots         = as.double(splineKnots),           
              m                   = as.integer(m - 1),
              bin                 = as.integer(0),
              v                   = as.double(0)          
    )  
    v <- out$v
  }
  return(v)
}

cubicISpline <- function(tau, interiorKnots, m, cIndicator = F){
  # Given knots splineKnots, evaluates the mth i-spline
  #  of 3 degrees of freedom at tau.
  # 
  # Args:
  #   tau: scalar in (0, 1).
  #   interiorKnots: interior knot sequence.
  #   m: index for the spline.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interiorKnots and splineDegreesOfFreedom degrees of freedom,
  #       the mth i-spline evaluated at tau.
  # Error handling
  checkUnitInterval(tau)
  checkInteriorKnots(interiorKnots)
  checkInteger(m, positive = T)
  splineKnots <- makeSplineKnots(interiorKnots, splineDegreesOfFreedom = 3, cIndicator = F)
  if(m  > length(splineKnots) - 3) {
    stop(message = "m cannot be greater than length(splineKnots) - 3.\n")
  }
  checkBinary(cIndicator)
  splineDegreesOfFreedom <- 3
  if(!cIndicator) {
    bin <- findInterval(tau,splineKnots)
    if(bin < m){
      v <- 0
    }
    else if ((bin - splineDegreesOfFreedom + 1) > m) {
      v <- 1
    }
    else if (bin == m){
      v <- ((tau - splineKnots[m]) ^ 3) /
        (
          (splineKnots[m+1]  - splineKnots[m]) *
            (splineKnots[m+2]  - splineKnots[m]) *
            (splineKnots[m+3]  - splineKnots[m])
        )
    }
    else if(bin == (m + 1)){
      i1 <- (tau - splineKnots[m]) / (splineKnots[m + 2] - splineKnots[m + 1])
      i2 <- ((tau - splineKnots[m + 1]) ^ 2 -
               (splineKnots[m + 3]-tau) ^ 2) /
        ((splineKnots[m + 3] - splineKnots[m + 1]) * (splineKnots[m + 2] - splineKnots[m + 1]))
      i3 <- (splineKnots[m + 3] - tau) ^ 3 /
        ((splineKnots[m + 3] - splineKnots[m + 1]) *
           (splineKnots[m + 3] - splineKnots[m]) *
           (splineKnots[m + 2]-splineKnots[m + 1])) -
        (tau - splineKnots[m]) ^ 3 /
        ((splineKnots[m + 3] - splineKnots[m]) *
           (splineKnots[m + 2] - splineKnots[m]) *
           (splineKnots[m+2]-splineKnots[m+1])
        )
      v <- i1 + i2 + i3
    }
    else{
      v <- 1 - (splineKnots[m + 3] - tau) ^ 3 /
        ((splineKnots[m + 3] - splineKnots[m + 2]) *
           (splineKnots[m + 3] - splineKnots[m + 1]) *
           (splineKnots[m + 3] - splineKnots[m]))
    }
  } else {
    out <- .C("cubicISpline",
              tau                 = as.double(tau),          
              splineKnotsLength   = as.integer(length(splineKnots)),          
              splineKnots         = as.double(splineKnots),          
              m                   = as.integer(m - 1),
              bin                 = as.integer(0),
              v                   = as.double(0)          
    )  
    v <- out$v
  }
  return(v)
}

makeCubicSplineCoefficients <- function(interiorKnots, cIndicator = F) {
  # Constructs array containing coefficients for cubic spline calculation.
  # Args:
  #   interiorKnots: a sequence of knots in (0, 1).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   splineCoefs: an array of dimension 4 x nBasisFunctions x (nBasisFunctions - 1)
  #     containing the coefficients for cubic spline calculation.
  #     splineCoefs[1, , ] contains the cubic coefficients.
  #     splineCoefs[2, , ] contains the quadratic coefficients.
  #     splineCoefs[3, , ] contains the linear coefficients.
  #     splineCoefs[4, , ] contains the constant coefficients.
  #     The second dimension indexes the basis function.
  #     The third dimension indexes the bin.
  #     Given that tau is in the mth bin,
  #       I_l(tau) =  splineCoefs[1,l, m] * tau ^ 3
  #           + splineCoefs[2,l, m] * tau ^ 2
  #           + splineCoefs[3,l, m] * tau
  #           + splineCoefs[4,l, m].
  # Note: This calculation begins with the equations in
  #   Zhou, Jingwen, Montserrat Fuentes, and Jerry Davis. "Calibration of
  #   numerical model output using nonparametric spatial density functions."
  #   Journal of agricultural, biological, and environmental statistics 16.4  
  #   (2011): 531-553.
  # a, b, and c are scaling factors
  # a represents the denominator in Equation 3.4,
  #   where gamma_m < tau < gamma_(m + 1).
  # b represents the denominators for i1*, i2*, and i3* (2 for i3*).
  # c represents the denominator in Equation 3.4,
  #   where gamma_(m + 2) < tau < gamma_(m + 3).
  # dA, dB, dC, and dD are (nBasisFunctions - 1) x (nBasisFunctions - 1) matrices.
  # row indexes basis function (spline, that is).
  # column indexes bin.
  # dA represents the coefficients for the cubic coefficient terms.
  # dB represents the quadratic coefficient.
  # dC represents the linear coefficient.
  # dD represents the constant coefficient.
  #Error Handling
  checkInteriorKnots(interiorKnots)
  checkBinary(cIndicator)
  splineKnots <- makeSplineKnots(interiorKnots)
  nBasisFunctions <- length(splineKnots) - 4
  if(!cIndicator) {
    a <- matrix(0, nrow = length(splineKnots), 1)
    b <- matrix(0, nrow = length(splineKnots), 4)
    c <- matrix(0, nrow = length(splineKnots), 1)
    for(m in 3:(nBasisFunctions - 1)) {
      a[m] <- 1 / ((splineKnots[m + 1] - splineKnots[m]) *
                     (splineKnots[m + 2] - splineKnots[m]) *
                     (splineKnots[m + 3] - splineKnots[m])
      )
    }
    for(bin in 3:(nBasisFunctions - 1)) {
      m <- bin - 1
      b[bin,1] <- 1 / (splineKnots[m + 2] - splineKnots[m + 1])
      b[bin,2] <- 1 /((splineKnots[m + 3] - splineKnots[m + 1])  * (splineKnots[m + 2] - splineKnots[m + 1]))
      b[bin,3] <- - 1/((splineKnots[m + 3] - splineKnots[m + 1]) * (splineKnots[m + 3] - splineKnots[m]) * (splineKnots[m + 2] - splineKnots[m + 1]))
      b[bin,4] <- - 1/((splineKnots[m + 3] - splineKnots[m]) * (splineKnots[m + 2] - splineKnots[m]) * (splineKnots[m + 2] - splineKnots[m + 1]))
    }
    for(bin in 3:(nBasisFunctions - 1)) {
      m <- bin - 2
      c[bin]   <-  1 / ((splineKnots[m + 3] - splineKnots[m + 2]) * (splineKnots[m + 3] - splineKnots[m + 1]) * (splineKnots[m+3] - splineKnots[m]))
    }
   
    dA <- dB <- dC <- dD <- matrix(0, nrow = nBasisFunctions - 1, ncol = nBasisFunctions - 1)
    for(l in 1:nBasisFunctions) { #indexes basis function
      for(m in 3: (nBasisFunctions - 1)) { #indexes bin
        if (m > (l + 2)) {dD[l,m] <- 1}
        if (l == m) {
          dA[l, m] <- a[m]
          dB[l, m] <- - 3 * a[m] * splineKnots[m]
          dC[l, m] <-   3 * a[m] * splineKnots[m] ^ 2
          dD[l, m] <- - a[m] * splineKnots[m] ^ 3
        }
        if ((l + 1) == m) {
          dA[l, m] <- b[m, 3] + b[m, 4]
          dB[l, m] <- - 3 * b[m, 3] * splineKnots[m + 2] - 3 * b[m, 4] * splineKnots[m - 1]
          dC[l, m] <-   3 * b[m, 3] * (splineKnots[m + 2] ^ 2) + 3 * b[m, 4] * (splineKnots[m - 1] ^ 2) + b[m,1] + 2 * b[m,2] * (splineKnots[m + 2] - splineKnots[m]) 
          dD[l, m] <- - b[m, 3] * splineKnots[m + 2] ^ 3 - b[m, 4] * (splineKnots[m - 1] ^ 3)  - b[m, 1] * splineKnots[m - 1]  + b[m, 2] * (splineKnots[m] ^ 2 - splineKnots[m + 2] ^ 2 )
        }
        if ((l + 2) == m) {
          dA[l, m] <- c[m]
          dB[l, m] <- -3 * c[m] * splineKnots[m + 1]
          dC[l, m] <- 3 * c[m] * splineKnots[m + 1]^2
          dD[l, m] <- -c[m] * splineKnots[m + 1] ^ 3 + 1
        }
      }
    }
   
    splineCoefs <- array(0, dim = c(4, nBasisFunctions, nBasisFunctions - 1))
    splineCoefs[1, 2 : nBasisFunctions,] <- dA
    splineCoefs[2, 2 : nBasisFunctions,] <- dB
    splineCoefs[3, 2 : nBasisFunctions,] <- dC
    splineCoefs[4, 2 : nBasisFunctions,] <- dD
   
    splineCoefs[4, 1, ] <- 1
  } else {
    out <- .C("makeCubicSplineCoefficients",
              nBasisFunctions = as.integer(nBasisFunctions),          
              splineKnots     = as.double(splineKnots),
              splineCoefs    = as.double(rep(0, 4 * nBasisFunctions * (nBasisFunctions - 1)))
    )     
    splineCoefs <- array(out$splineCoefs, dim = c(4, nBasisFunctions, nBasisFunctions - 1))
  }
  return(splineCoefs)
}

cubicMSpline2 <- function(tau, splineKnots, splineCoefs, m,  cIndicator = F) {
  # Given knots splineKnots and splineCoefs, evaluates the mth m-spline
  #  of 3 degrees of freedom at tau.
  # 
  # Args:
  #   tau: scalar in (0, 1).
  #   splineKnots: knot sequence.
  #   splineCoefs: coefficients used to calculate the cubic spline.
  #      Assumed to be correctly calculated inside of this function.
  #   m: index for the spline.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interiorKnots and splineDegreesOfFreedom degrees of freedom,
  #       the mth m-spline evaluated at tau.
  # Error handling
  checkUnitInterval(tau)
  checkCubicSplineKnots(splineKnots)
  nBasisFunctions <- length(splineKnots) - 4
  checkNumeric(splineCoefs, lengthOne = F)
  if(any(dim(splineCoefs) != c(4, nBasisFunctions, nBasisFunctions - 1))) {
    stop(message = "splineCoefs must be a 4 x nBasisFunctions x nBasisFunctions - 1 array.\n")
  }
  checkInteger(m, positive = T)
  if(m  > length(splineKnots) - 3) {
    stop(message = "m cannot be greater than length(splineKnots) - 3.\n")
  }
  checkBinary(cIndicator)
  splineDegreesOfFreedom <- 3
  if(!cIndicator) {
    bin <- findInterval(tau, splineKnots)
    a <- splineCoefs[1, m + 1, bin]
    b <- splineCoefs[2, m + 1, bin]
    c <- splineCoefs[3, m + 1, bin]
    d <- splineCoefs[4, m + 1, bin]
    v <- 3 * a * tau ^ 2 + 2 * b * tau  + c
  } else {
    out <- .C("cubicMSpline2",
              tau               = as.double(tau),
              splineKnotsLength = as.integer(length(splineKnots)),          
              splineKnots       = as.double(splineKnots),
              splineCoefs       = as.double(splineCoefs),
              m                 = as.integer(m - 1),
              bin               = as.integer(2),
              v                 = as.double(0)
    )     
    v <- out$v
  }
  return(v)
}  

cubicISpline2 <- function(tau, splineKnots, splineCoefs, m, cIndicator = F) {
  # Given knots splineKnots and splineCoefs, evaluates the mth i-spline
  #  of 3 degrees of freedom at tau.
  # 
  # Args:
  #   tau: scalar in (0, 1).
  #   splineKnots: knot sequence.
  #   splineCoefs: coefficients used to calculate the cubic spline.
  #      Assumed to be correctly calculated inside of this function.
  #   m: index for the spline.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interiorKnots and splineDegreesOfFreedom degrees of freedom,
  #       the mth i-spline evaluated at tau.
  # Error handling:
  checkUnitInterval(tau)
  checkCubicSplineKnots(splineKnots)
  nBasisFunctions <- length(splineKnots) - 4
  checkNumeric(splineCoefs, lengthOne = F)
  if(any(dim(splineCoefs) != c(4, nBasisFunctions, nBasisFunctions - 1))) {
    stop(message = "splineCoefs must be a 4 x nBasisFunctions x nBasisFunctions - 1 array.\n")
  }
  checkInteger(m, positive = T)
  if(m  > length(splineKnots) - 3) {
    stop(message = "m cannot be greater than length(splineKnots) - 3.\n")
  }
  checkBinary(cIndicator)
  splineDegreesOfFreedom <- 3
  if(!cIndicator) {
    bin <- findInterval(tau, splineKnots)
    a <- splineCoefs[1, m + 1, bin]
    b <- splineCoefs[2, m + 1, bin]
    c <- splineCoefs[3, m + 1, bin]
    d <- splineCoefs[4, m + 1, bin]
    v <- a * tau ^ 3 + b * tau ^ 2 + c * tau + d
  } else {
    out <- .C("cubicISpline2",
              tau               = as.double(tau),
              splineKnotsLength = as.integer(length(splineKnots)),          
              splineKnots       = as.double(splineKnots),
              splineCoefs       = as.double(splineCoefs),
              m                 = as.integer(m - 1),
              bin               = as.integer(2),
              v                 = as.double(0)
    )     
    v <- out$v
  }
  return(v)
}  

makeStickBreakingWeights <- function(v, cIndicator = F){
  # Transforms a vector of beta random variables and a 1 at the end into
  #  latent class weights.
  # 
  # Args:
  #   v: vector containing beta random variables,
  #      except for the final element which is 1.
  #    
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   nu : Latent class weights.
  # Error handling:
  checkNumeric(v, lengthOne = F)
  lengthV <- length(v)
  if(v[lengthV] != 1) {
    stop(message = "The last element of v must be 1.\n")
  }
  if(lengthV > 1) {
    checkUnitInterval(v[1 : (lengthV - 1)], lengthOne = F)
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    if(lengthV == 1) {
      nu <- v
    } else {
      nu <- v
      nu[2 : lengthV] <- nu[2 : lengthV] * cumprod(1 - v[1:(lengthV - 1)])
    }
  } else {
    out <- .C("makeStickBreakingWeights",
              lengthV  = as.integer(lengthV),
              v         = as.double(v),
              nu        = as.double(rep(0, lengthV))          
    )     
    nu <- out$nu
  }
  return(nu)
}

logDensityNormal <- function(y, loc, scale, shape = 0, cIndicator = F) {
  # Finds fY(y), where Y is normally distributed with
  ## location scale parameters loc/scale
  #
  # Args:
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter (ignored).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    densityY <- dnorm(y, loc, sd = scale, log = T)
  } else {
    out <- .C("logDensityNormal",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              densityY = as.double(0)
    )  
    densityY <- out$densityY
  }
  return(densityY)
}

logDensityT <- function(y, loc, scale, shape, cIndicator = F) {
  # Finds fY(y), where Y is student's t distributed with
  ## location/scale/df parameters loc/scale/shape
  #
  # Args:
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    z <- (y - loc) / scale
    densityY <- dt(z, df = shape, log = T) - log(scale)
  } else {
    out <- .C("logDensityT",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              densityY = as.double(0)
    )
    densityY <- out$densityY
  }
  return(densityY)
}

logDLogistic <- function(y, loc, scale, shape = 0, cIndicator = F) {
  # Finds fY(y), where Y is logistically distributed with
  ## location scale parameters loc/scale
  #
  # Args:
  #   y: scalar
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter (ignored).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    densityY <- dlogis(y, loc, scale, log = T)
  } else {
    out <- .C("logDLogistic",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              densityY = as.double(0)
    )
    densityY <- out$densityY  
  }
  return(densityY)
}

logDensityAsymmetricLaplace <- function(y, loc, scale, shape, cIndicator = F) {
  # Finds fY(y), where Y is asymmetric Lapace distributed with
  ## location/scale/shape parameters loc/scale/shape
  #
  # Args:
  #   cIndicator: an indicator if the computation should be performed in C.
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter in (0, 1).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    densityY <- log(dalap(y, loc, scale = scale, tau = shape))
  } else {
    out <- .C("logDensityAsymmetricLaplace",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              densityY = as.double(0)
    )
    densityY <- out$densityY
  }
  return(densityY)
}

logDensityWeibull <- function(y, loc = 0, scale, shape, cIndicator = F) {
  # Finds fY(y), where Y is Weibull distributed with
  ## scale/shape parameters scale/shape.
  #
  # Args:
  #   y: scalar.
  #   loc: location parameter (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  checkNumeric(y)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    densityY <- dweibull(y, shape, scale, log = T)
  } else {
    out <- .C("logDensityWeibull",
              y = as.double(y),
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              densityY = as.double(0)
    )
    densityY <- out$densityY
  }
  return(densityY)
}

logDensityGamma <- function(y, loc = 0, scale, shape, cIndicator = F) {
  # Finds fY(y), where Y is gamma distributed with
  ## scale/shape parameters scale/shape.
  #
  # Args:
  #   y: scalar.
  #   loc: positive scale parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  checkNumeric(y)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    densityY <- dgamma(y, shape = shape, scale = scale, log = T)
  } else {
    out <- .C("logDensityGamma",
              y = as.double(y),
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              densityY = as.double(0)
    )
    densityY <- out$densityY
  }
  return(densityY)
}

logDensityBeta <- function(y, shape1, shape2, cIndicator = F) {
  # Finds fY(y), where Y is beta distributed with
  ## parameters shape1 and shape2.
  #
  # Args:
  #   y: scalar.
  #   shape1: positive first shape parameter.
  #   shape1: positive second shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  checkNumeric(y)
  checkNumeric(shape1, positive = T)
  checkNumeric(shape2, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    densityY <- dbeta(y, shape1, shape2, log = T)
  } else {
    out <- .C("logDensityBeta",
              y = as.double(y),
              shape1  = as.double(shape1),
              shape2  = as.double(shape2),
              densityY = as.double(0)
    )
    densityY <- out$densityY
  }
  return(densityY)
}

logDensityInverseWishart <- function(df, sai, x, cIndicator = F) {
  # evaluates the log density of the inverse wishart distribution evaluated at x
  #   with scale matrix sai and df degrees of freedom.
  # Args:
  #   df = degrees of freedom.
  #   sai = scale matrix.
  #   x = positive definite matrix argument.
  # Returns:
  #  densityY: log density.
  # Error handling:
  checkNumeric(df, positive = T)
  checkInvertibleSymmetricMatrix(sai)
  checkInvertibleSymmetricMatrix(x)
  p <- nrow(x)
  if(nrow(sai) != p) {
    stop(message = "sai must have the same number of rows as x.\n")
  }
  if(df <= (p - 1)) {
    stop(message = "df must be greater than p - 1.\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    densityY <- log(diwish(x, df, sai))
  } else {
    out <- .C("logDensityInverseWishart",
              p = as.integer(p),
              df = as.double(df),
              sai = as.double(sai),
              x = as.double(x),
              densityY = as.double(0)
    )
    densityY <- out$densityY
  }
  return(densityY)
}  

cdfNormal <- function(y, loc, scale, shape = 0, cIndicator = F) {
  # Finds P(Y <= y), where Y is normally distributed with
  ## location scale parameters loc/scale.
  #
  # Args:
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    cdfY <- pnorm(y, loc, sd = scale)
  } else {
    out <- .C("cdfNormal",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              cdfY = as.double(0)
    )
    cdfY <- out$cdfY
  }
  return(cdfY)
}

cdfT <- function(y, loc, scale, shape, cIndicator = F) {
  # Finds P(Y <= y), where Y is student's t distributed with
  ## location/scale/df parameters loc/scale/shape.
  #
  # Args:
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    cdfY <- pt((y - loc) / scale, df = shape)
  } else {
    out <- .C("cdfT",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              cdfY = as.double(0)
    )
    cdfY <- out$cdfY
  }
  return(cdfY)
}

cdfLogistic <- function(y, loc, scale, shape = 0, cIndicator = F) {
  # Finds P(Y <= y), where Y is logistically distributed with
  ## location scale parameters loc/scale.
  #
  # Args:
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    cdfY <- plogis(y, loc, scale)
  } else {
    out <- .C("cdfLogistic",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              cdfY = as.double(0)
    )
    cdfY <- out$cdfY
  }
  return(cdfY)
}

cdfAsymmetricLaplace <- function(y, loc, scale, shape, cIndicator = F) {
  # Finds P(Y <= y), where Y is asymmetric Lapace distributed with
  ## location/scale/shape parameters loc/scale/shape.
  #
  # Args:
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter in (0, 1).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  checkNumeric(y)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    cdfY <- palap(y, loc, scale, tau = shape)
  } else {
    out <- .C("cdfAsymmetricLaplace",
              y = as.double(y),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              cdfY = as.double(0)
    )
    cdfY <- out$cdfY
  }
  return(cdfY)
}

cdfWeibull <- function(y, loc = 0, scale, shape, cIndicator = F) {
  # Finds P(Y <= y), where Y is Weibull distributed with
  ## scale/shape parameters scale/shape,
  #
  # Args:
  #   y: scalar.
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y)
  # Error handling:
  checkNumeric(y)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    cdfY <- pweibull(y, shape, scale)
  } else {
    out <- .C("cdfWeibull",
              y = as.double(y),
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              cdfY = as.double(0)
    )
    cdfY <- out$cdfY
  }
  return(cdfY)
}

pGamma <- function(y, loc, scale, shape, cIndicator = F) {
  # Finds P(Y <= y), where Y is gamma distributed with
  ## scale/shape parameters scale/shape.
  #
  # Args:
  #   y: scalar.
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  checkNumeric(y)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    cdfY <- pgamma(y, shape = shape, scale = scale)
  } else {
    out <- .C("pGamma",
              y = as.double(y),
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              cdfY = as.double(0)
    )
    cdfY <- out$cdfY
  }
  return(cdfY)
}

qNorm <- function(tau, loc, scale, shape = 0, cIndicator = F) {
  # Finds the tauth quantile of Y, where Y is normally distributed with
  ## location scale parameters loc/scale.
  #
  # Args:
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value qTau, where P(Y <= qTau) = tau.
  # Error handling:
  checkUnitInterval(tau)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    qTau <- qnorm(tau, loc, sd = scale)
  } else {
    out <- .C("qNorm",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              qTau = as.double(0)
    )
    qTau <- out$qTau
  }
  return(qTau)
}

qT <- function(tau, loc, scale, shape, cIndicator = F) {
  # Finds the tauth quantile of Y, where Y is student's t distributed with
  ## location/scale/df parameters loc/scale/shape.
  #
  # Args:
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value qTau, where P(Y <= qTau) = tau.
  # Error handling:
  checkUnitInterval(tau)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    qTau <- scale * qt(tau, df = shape) + loc
  } else {
    out <- .C("qT",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              qTau = as.double(0)
    )
    qTau <- out$qTau
  }
  return(qTau)
}

qLogistic <- function(tau, loc, scale, shape = 0, cIndicator = F) {
  # Finds the tauth quantile of Y, where Y is logistically distributed with
  ## location scale parameters loc/scale.
  #
  # Args:
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value qTau, where P(Y <= qTau) = tau.
  # Error handling:
  checkUnitInterval(tau, lengthOne = F)
  checkNumeric(loc, lengthOne = F)
  checkNumeric(scale, lengthOne = F, positive = T)
  n <- length(tau)
  if(length(loc) == 1) {
    loc <- rep(loc, n)
  }
  if(length(scale) == 1) {
    scale <- rep(scale, n)
  }
  if(length(loc) != n) {
    stop(message = "tau and loc must be of the same length.\n")
  }
  if(length(scale) != n) {
    stop(message = "tau and scale must be of the same length.\n")
  }
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    qTau <- qlogis(tau, loc, scale)
  } else {
    if(n != 1) {
      stop(message = "C version of this function assumes scalar inputs.\n")
    }
    out <- .C("qLogistic",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              qTau = as.double(0)
    )
    qTau <- out$qTau
  }
  return(qTau)
}

qAsymmetricLaplace <- function(tau, loc, scale, shape, cIndicator = F) {
  # Finds the tauth quantile of Y, where Y is asymmetric Lapace distributed with
  ## location/scale/shape parameters loc/scale/shape.
  #
  # Args:
  #   cIndicator: an indicator if the computation should be performed in C.
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter in (0, 1).
  # Returns:
  #   The value qTau, where P(Y <= qTau) = tau.
  # Error handling:
  checkUnitInterval(tau)
  checkNumeric(loc)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    qTau <- qalap(tau, loc, scale, tau = shape)
  } else {
    out <- .C("qAsymmetricLaplace",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              qTau = as.double(0)
    )
    qTau <- out$qTau
  }
  return(qTau)
}

qWeibull <- function(tau, loc, scale, shape, cIndicator = F) {
  # Finds the tauth quantile of Y, where Y is Weibull distributed with
  ## scale/shape parameters scale/shape.
  #
  # Args:
  #   tau: quantile level in (0, 1).
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value qTau, where P(Y <= qTau) = tau.
  # Error handling:
  checkUnitInterval(tau)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    qTau <- qweibull(tau, shape, scale)
  } else {
    out <- .C("qWeibull",
              tau = as.double(tau),
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              qTau = as.double(0)
    )
    qTau <- out$qTau
  }
  return(qTau)
}

qGamma <- function(tau, loc, scale, shape, cIndicator = F) {
  # Finds the tauth quantile of Y, where Y is gamma distributed with
  ## scale/shape parameters scale/shape.
  #
  # Args:
  #   tau: quantile level in (0, 1).
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value qTau, where P(Y <= qTau) = tau.
  # Error handling:
  checkUnitInterval(tau)
  checkNumeric(scale, positive = T)
  checkNumeric(shape, positive = T)
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    qTau <- qgamma(tau, shape = shape, scale = scale)
  } else {
    out <- .C("qGamma",
              tau = as.double(tau),
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              qTau = as.double(0)
    )
    qTau <- out$qTau
  }
  return(qTau)
}

truncatedGaussianMoments <- function(mu, sigma, a, b, type, cIndicator = F) {
  # Finds mean and variance for a truncated Gaussian random variable.
  # Args:
  #   mu: mean.
  #   sigma: standard deviation.
  #   a: left-censoring point (can be -Inf).
  #   b: right-censoring point (can be Inf).
  #   type: an indicator for the type of censoring.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   mn:  means of the distribution.
  #   vr:  variances of the distribution.
  # Error handling
 
  checkNumeric(mu)
  checkNumeric(sigma, positive = T)
  checkNumeric(a)
  checkNumeric(b)
  if(b <= a) {
    stop(message = "b must be greater than a.\n")
  }
 
  if(length(type) != 1) {
    stop(message = "type must be of length 1.\n")
  }
  typecheck <- type %in% c("left", "interval", "right")
  if(!typecheck) {
    stop(message = "type must be equal to \"left\", \"interval\" or \"right\".\n")
  }
 
  if(type == "left" & a == -Inf) {
    stop(message = "a cannot be -Inf if type = left.\n")
  }
 
  if(type == "right" & b == Inf) {
    stop(message = "b cannot be Inf if type = right.\n")
  }
 
  if(type == "interval" & a == -Inf) {
    stop(message = "a cannot be -Inf if type = interval.\n")
  }
  if(type == "interval" & b == Inf) {
    stop(message = "b cannot be Inf if type = interval.\n")
  }
 
  sigma2 <- sigma * sigma
 
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    if(type == "left") {
      aCentered <- (a - mu) / sigma
      lambda <- dnorm(aCentered) / (1 - pnorm(aCentered))
      delta <- lambda * (lambda - aCentered)
      mn <- mu + sigma * lambda
      vr <- sigma2 * (1 - delta)    
    } else if (type == "right") {
      bCentered <- (b - mu) / sigma
      bCenteredRatio <- dnorm(bCentered) / pnorm(bCentered)
      mn <- mu - sigma * bCenteredRatio
      vr <- sigma2 * (1 - bCentered * bCenteredRatio - bCenteredRatio ^ 2)    
    } else {
      aCentered <- (a - mu) / sigma
      bCentered <- (b - mu) / sigma  
      ratio <- (dnorm(aCentered) - dnorm(bCentered)) / (pnorm(bCentered) - pnorm(aCentered))    
      mn <- mu + sigma * ratio
      vr <- sigma2 * (1 +
                        (aCentered * dnorm(aCentered) - bCentered * dnorm(bCentered)) /
                        (pnorm(bCentered) - pnorm(aCentered)) -
                        ratio ^ 2    
      )
    }
  } else {
    if(type == "left") {
      type <- 1
    }
    if(type == "right") {
      type <- 2
    }
    if(type == "interval") {
      type <- 3
    }
    a <- ifelse(a == -Inf, -99, a)
    b <- ifelse(b ==  Inf, -99, b)  
    out <- .C("truncatedGaussianMoments",
              type        = as.integer(type),
              mu          = as.double(mu),
              sigma       = as.double(sigma),
              a           = as.double(a),        
              b           = as.double(b),
              mn          = as.double(0),
              vr          = as.double(0)
    )          
    mn <- out$mn
    vr <- out$vr
  }
  list(mn = mn, v = vr)
}

randomNormal <- function(nSamples, mu, sigma, cIndicator = F) {
  # Generates a multivariate normal random variable.
  # The sample will differ depending on cIndicator, but both versions are valid.
  #
  # Args:
  #   nSamples: number of samples.
  #   mu: the mean parameter.
  #   sigma: the covariance parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   y: multivariate normal realization.
  #
  # Error handling:
  checkInteger(nSamples, positive = T)
  checkNumeric(mu, lengthOne = F)
  checkInvertibleSymmetricMatrix(sigma)
  dimY <- length(mu)
  if(nrow(sigma) != dimY) {
    stop(message = "number of rows of sigma and length of mu must be equal.\n")
  }
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    y <- rmvnorm(nSamples, mean = mu, sigma = sigma)
  } else {
    out <- .C("randomNormal",
              nSamples   = as.integer(nSamples),
              dimY       = as.integer(dimY),          
              mu          = as.double(mu),          
              sigma       = as.double(sigma),          
              y           = as.double(rep(0,  nSamples * dimY))
    )  
    y <- array(out$y, dim = c(nSamples, dimY))
  }
  return(y)
}

randomTruncatedNormal <- function(nSamples, mu, sigma, a, b, cIndicator = F){
  # Generates a truncated normal random variable.
  #
  # Args:
  #   nSamples: number of samples.
  #   mu: the mean parameter.
  #   sigma: the scale parameter.
  #   a: the left truncation point.  Can be -Inf.
  #   b: the right truncation point.  Can be Inf.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   y: truncated normal realization.
  #
  # Error handling:
  checkInteger(nSamples, positive = T)
  checkNumeric(mu, lengthOne = F)
  if(length(mu) == 1) {
    mu <- rep(mu, nSamples)
  }
  if(length(mu) != nSamples) {
    stop(message = "length of mu must be nSamples.\n")
  }
 
  checkNumeric(sigma, positive = T, lengthOne = F)
  if(length(sigma) == 1) {
    sigma <- rep(sigma, nSamples)
  }
  if(length(sigma) != nSamples) {
    stop(message = "length of sigma must be nSamples.\n")
  }
  checkNumeric(a, lengthOne = F)
  if(length(a) == 1) {
    a <- rep(a, nSamples)
  }
  if(length(a) != nSamples) {
    stop(message = "length of a must be nSamples.\n")
  }
  checkNumeric(b, lengthOne = F)
  if(length(b) == 1) {
    b <- rep(b, nSamples)
  }
  if(length(b) != nSamples) {
    stop(message = "length of b must be nSamples.\n")
  }
  if(any(b < a)) {
    stop(message = "b must be greater than or equal to a for all elements.\n")
  }
  checkBinary(cIndicator)
 
  if(!cIndicator) {
    u <- runif(nSamples)
    pNormAlpha <- ifelse(a == -Inf, 0, pnorm((a - mu) / sigma))
    pNormBeta <- ifelse(b == Inf, 1, pnorm((b - mu) / sigma))
    y <- qnorm(pNormAlpha + u * (pNormBeta - pNormAlpha)) * sigma + mu
  } else {
    a <- ifelse(a == -Inf, -99, a)
    b <- ifelse(b ==  Inf, -99, b)  
    out <- .C("randomTruncatedNormal",
              nSamples    = as.integer(nSamples),
              mu          = as.double(mu),
              sigma       = as.double(sigma),
              a           = as.double(a),        
              b           = as.double(b),
              y           = as.double(rep(0, nSamples))
    )          
    y <- matrix(out$y, ncol = 1)
  }
  return(y)
}

randomInverseWishart <- function(df, sai, cIndicator = F) {
  # Generates an inverse Wishart random variable.
  # The sample will differ depending on cIndicator, but both versions are valid.
  #
  # Args:
  #   df: the degrees of freedom.
  #   sai: the scale matrix.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   w: inverse Wishart realization.
  #
  # Error handling:
  checkNumeric(df, positive = T)
  checkInvertibleSymmetricMatrix(sai)
  checkBinary(cIndicator)
  if(!cIndicator) {
    w <- riwish(df, sai)
  } else {
    nRow <- nrow(sai)
    out <- .C("randomInverseWishart",
              nRow    = as.integer(nRow),
              df      = as.double(df),        
              sai     = as.double(sai),
              w       = as.double(rep(0, nRow * nRow))
    )          
    w <- array(out$w, dim = c(nRow, nRow))
  }
  return(w)
}

randomDirichlet <- function(nSamples = 1, probs, cIndicator = F) {
  # Generates nSamples draws from a Dirichlet(probs) distribution.
  # 
  # Args:
  #   nSamples: the number of samples.
  #   probs: vector of probabilities.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : nSamples x length(probs) array of Dirichlet draws.
  # Error handling:
  checkInteger(nSamples, positive = T)
  checkUnitInterval(probs, lengthOne = F)
  checkBinary(cIndicator)
  if(!cIndicator) {
    x <- rdirichlet(nSamples, probs)
  } else {
    nCells <- length(probs)  
    out <- .C("randomDirichlet",
              nSamples       = as.integer(nSamples),          
              nCells         = as.integer(nCells),          
              probs          = as.double(probs),
              x              = as.double(rep(0, nSamples * nCells))
    )
    x <- array(out$x, dim = c(nSamples, nCells))
  }
  return(x)
}

randomMultinomial <- function(nSamples = 1, probs, cIndicator = F) {
  # Generates nSamples draws from a multinomial(probs) distribution.
  # 
  # Args:
  #   nSamples: the number of samples.
  #   probs: vector of probabilities.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : length(probs) vector of multinomial draws.
  # Error handling:
  checkInteger(nSamples, positive = T)
  checkUnitInterval(probs, lengthOne = F)
  checkBinary(cIndicator)
  if(!cIndicator) {
    x <- rmultinom(n = 1, size = nSamples, prob = probs)
  } else {
    nCells <- length(probs)  
    out <- .C("randomMultinomial",
              nSamples         = as.integer(nSamples),          
              nCells           = as.integer(nCells),
              probs             = as.double(probs),
              x                 = as.integer(rep(0, nCells))          
    )
    x <- matrix(out$x, ncol = 1)
  }
  return(x)
}

randomUniform <- function(a = 0, b = 1, cIndicator = F) {
  # Generates a draw from a unif(a, b) distribution.
  # 
  # Args:
  #   a: lower limit in the support.
  #   b: upper limit in the support.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : unif(a, b) variable.
  # Error handling:
  checkNumeric(a)
  checkNumeric(b)
  if(a > b) {
    stop(message = "a must be less than or equal to b.\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    x <- runif(1, a, b)
  } else {
    out <- .C("randomUniform",
              a = as.double(a),
              b = as.double(b),
              x = as.double(0)          
    )
    x <- out$x
  }
  return(x)
}

randomBeta <- function(shape1, shape2, cIndicator = F) {
  # Generates a draw from a beta(shape1, shape2) distribution.
  # 
  # Args:
  #   shape1: the first shape parameter.
  #   shape2: the second shape parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : beta random variable.
  # Error handling:
  checkNumeric(shape1, positive = T)
  checkNumeric(shape2, positive = T)
  checkBinary(cIndicator)
  if(!cIndicator) {
    x <- rbeta(n = 1, shape1, shape2)
  } else {
    out <- .C("randomBeta",
              shape1 = as.double(shape1),
              shape2 = as.double(shape2),
              x      = as.double(0)          
    )
    x <- out$x
  }
  return(x)
}

randomGamma <- function(shape, scale, cIndicator = F) {
  # Generates a draw from a gamma(shape, scale) distribution.
  # 
  # Args:
  #   shape: the shape parameter.
  #   scale: the scale parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : gamma random variable.
  # Error handling:
  checkNumeric(shape, positive = T)
  checkNumeric(scale, positive = T)
  checkBinary(cIndicator)
  if(!cIndicator) {
    x <- rgamma(n = 1, shape = shape, scale = scale)
  } else {
    nCells <- length(probs)  
    out <- .C("randomGamma",
              shape = as.double(shape),
              scale = as.double(scale),
              x     = as.double(0)          
    )
    x <- out$x
  }
  return(x)
}

randomExponential <- function(scale, cIndicator = F) {
  # Generates a draw from a gamma(shape, scale) distribution.
  # 
  # Args:
  #   scale: the scale parameter.
  #   cIndicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : exponential random variable.
  # Error handling:
  checkNumeric(scale, positive = T)
  checkBinary(cIndicator)
  if(!cIndicator) {
    x <- rexp(n = 1, rate = 1 / scale)
  } else {
    out <- .C("randomExponential",
              scale = as.double(scale),
              x     = as.double(0)          
    )
    x <- out$x
  }
  return(x)
}
  
updateKroneckerInverseWishart <- function(nu0, sigma0, e, omega, cIndicator = F) {
  # Performs a Gibbs update of the posterior inverse Wishart random variable sigma,
  # where e|omega, sigma ~ N(0, calculateKronecker(omega, sigma)) and sigma ~ IW(nu0, sigma0).
  #
  # Args:interiorKnots_length
  #   nRowOmega: the number of rows of omega.
  #   nRowSigma: the number of rows of sigma.
  #   nu0: the prior scale for sigma.
  #   sigma0: the prior location for sigma.
  #   e: a multivariate normal random variable that has mean 0
  #     and covariance calculateKronecker(omega ^ -1, sigma)).
  #   omega: a valid precision matrix.
  #
  # Returns:
  #   sigma: the updated covariance.
  # Error Handling:
  checkNumeric(nu0, positive = T)
  checkInvertibleSymmetricMatrix(sigma0)
  checkNumeric(e, lengthOne = F)
  if(!isSymmetric(omega)) {
    stop(message = "omega should be a symmetric matrix.\n")
  }
 
  nRowOmega <- nrow(omega)
  nRowSigma <- nrow(sigma0)
 
  if(length(e) != (nRowOmega * nRowSigma)) {
    stop(message = "e must be of length nrow(omega) * nrow(sigma0).\n")
  }
  checkBinary(cIndicator)
  if(!cIndicator) {
    s <- sigma0
    for(i in 1:nRowOmega) {
      for(j in 1:nRowOmega) {
        thisE1 <- e[((i - 1) * nRowSigma + 1) : (i * nRowSigma)]
        thisE2 <- e[((j - 1) * nRowSigma + 1) : (j * nRowSigma)]
        s <- s + omega[i,j] * (thisE2 %*% t(thisE1))
      }
    }
    sigma <- riwish(nRowOmega + nu0, s)
  } else {
    out <- .C("updateKroneckerInverseWishart",
              nRowOmega      = as.integer(nRowOmega),
              nRowSigma      = as.integer(nRowSigma),
              nu0            = as.double(nu0),
              sigma0         = as.double(sigma0),
              e              = as.double(e),
              omega          = as.double(omega),
              sigma          = as.double(rep(0, nRowSigma * nRowSigma))
             
    )
    sigma <- array(out$sigma, dim = c(nRowSigma, nRowSigma))
  }
  return(sigma)
}  
  
  
