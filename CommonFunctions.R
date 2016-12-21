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
## 1) CheckCharacter
## 2) CheckBinary
## 3) CheckNumericValues
## 4) CheckScalar
## 6) CheckNonnegativeNumber
## 7) CheckPositiveNumber
## 8) CheckPositiveScalar
## 9) CheckInteger
## 10) CheckUnitInverval
## 11) CheckLocs
## 12) CheckInteriorKnots
## 13) CheckInvertibileSymmetricMatrix

## 14) SetEqual
## 14) Mean
## 15) Variance
## 16) FindIntegerMax
## 16) FindInterval
## 17) MatrixMultiply
## 18) Kronecker
## 19) QuadraticForm
## 20) QuadraticFormKronecker
## 21) CholMatrix
## 22) InvertSymmetricMatrix
## 23) MakeAbsoluteDistance
## 24) MakeAutoregressivePrecision
## 25) MakeSpatialDistance
## 26) MakeSpatialPrecision

## 27) MakeSplineKnots
## 28) MSpline
## 29) ISpline
## 30) CubicMSpline
## 31) CubicISpline
## 32) MakeCubicSplineCoefficients
## 33) CubicMSpline2
## 34) CubicISpline2
## 35) MakeStickBreakingWeights

## 36) LogDNorm
## 37) LogDT
## 38) LogDLogistic
## 39) LogDAsymmetricLapace
## 40) LogDWeibull
## 41) LogDGamma
## 42) PNorm
## 43) PT
## 44) PLogistic
## 45) PAsymmetricLapace
## 46) PWeibull
## 47) PGamma
## 48) QNorm
## 49) QT
## 50) QLogistic
## 51) QAsymmetricLapace
## 52) QWeibull
## 53) QGamma
## 54) TruncatedGaussianMoments

## 55) RandomNormal
## 56) RandomTruncatedNormal
## 57) RandomInverseWishart
## 58) RandomDirichlet
## 59) RandomMultinomial
## 60) RandomBeta
## 61) RandomGamma
## 62) UpdateKroneckerInverseWishart

CheckCharacter <- function(v, v.name = NULL) {
  # Checks that an argument is character.
  #
  # Args:
  #   v: an object to be tested if binary.
  #   v.name: the name of v passed along from another function. 
  # Returns nothing.
  #
  # Error handling:
  v.name <- ifelse(is.null(v.name), deparse(substitute(v)), v.name)
  if(!is.character(v)) {
    stop(message = 
           paste0(v.name, " must be a character.\n"))
  }  
}

CheckBinary <- function(v, v.name = NULL, no.missing = T, length.one = T) {
  # Checks that an argument is binary 0/1.
  #
  # Args:
  #   v: an object to be tested if binary.
  #   v.name: the name of v passed along from another function. 
  #   no.missing: an indicator that no missing values are allowed.
  #   length.one: an indicator that the argument is of length 1.
  # Returns nothing.
  #
  # Error handling:
  if(length(no.missing) != 1) {
    stop(message = "no.missing must be of length 1.\n")
  }
  check1 <-  no.missing %in% c(T, F)
  if(!check1) {
    stop(message = "no.missing must be either TRUE or FALSE.\n")
  }
  if(length(length.one) != 1) {
    stop(message = "length.one must be of length 1.\n")
  }
  check2 <-  length.one %in% c(T, F)
  if(!check2) {
    stop(message = "length.one must be either TRUE or FALSE.\n")
  }
  v.name <- ifelse(is.null(v.name), deparse(substitute(v)), v.name)
  CheckCharacter(v.name)
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(v.name, " cannot contain only NA.\n"))
  }
  if(no.missing) {
    if(any(is.na(v))) {
      stop(message = 
             paste0(v.name, " must not have missing values.\n"))
    }
  }
  if(length.one) {
    if(length(v) != 1) {
      stop(message = paste0(v.name, " must be of length 1.\n"))
    }    
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
  zero.test <- v == 0 
  one.test <- v == 1
  if(any(!zero.test & !one.test)) {
    stop(message = paste0(v.name, " must be composed of 0 or 1.\n"))
  }
}

CheckInteger <- function(v, v.name = NULL, no.missing = T, length.one = T, 
                         nonnegative = F, positive = F) {
  # Checks that a scalar argument is an integer.
  #  It does not have to be of integer type.  
  # Args: 
  #   v: an object to be tested if numeric.
  #   v.name: the name of v passed along from another function. 
  #   no.missing: an indicator that no missing values are allowed.
  #   length.one: an indicator that the argument is of length 1.
  #   nonnegative: an indicator that values should be nonnegative.
  #   positive: an indicator that values should be positive.
  # Returns nothing.
  #
  # Error handling:
  v.name <- ifelse(is.null(v.name), deparse(substitute(v)), v.name)
  CheckCharacter(v.name)
  CheckBinary(no.missing)
  CheckBinary(length.one)
  CheckBinary(nonnegative)
  CheckBinary(positive)
  
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(v.name, " cannot contain only NA.\n"))
  }
  if(no.missing) {
    if(any(is.na(v))) {
      stop(message = 
             paste0(v.name, " must not have missing values.\n"))
    }
  }
  if(length.one) {
    if(length(v) != 1) {
      stop(message = paste0(v.name, " must be of length 1.\n"))
    }    
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
  if(any(v %% 1 != 0)) {
    stop(message = 
           paste0(v.name, " must be an integer.\n"))
  }
  if(nonnegative) {
    if(any(v < 0)) {
      stop(message = 
             paste0(v.name, " must have nonnegative elements.\n"))
    }
  }
  if(positive) {
    if(any(v <= 0)) {
      stop(message = 
             paste0(v.name, " must have positive elements.\n"))
    }
  }
}

CheckNumeric <- function(v, v.name = NULL, no.missing = T, length.one = T, 
                               nonnegative = F, positive = F) {
  # Checks that an argument is numeric.
  #
  # Args:
  #   v: an object to be tested if numeric.
  #   v.name: the name of v passed along from another function. 
  #   no.missing: an indicator that no missing values are allowed.
  #   length.one: an indicator that the argument is of length 1.
  #   nonnegative: an indicator that values should be nonnegative.
  #   positive: an indicator that values should be positive.
  # Returns nothing.
  #
  # Error handling:
  v.name <- ifelse(is.null(v.name), deparse(substitute(v)), v.name)
  CheckCharacter(v.name)
  CheckBinary(no.missing)
  CheckBinary(length.one)
  CheckBinary(nonnegative)
  CheckBinary(positive)
  
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(v.name, " cannot contain only NA.\n"))
  }
  if(no.missing) {
    if(any(is.na(v))) {
      stop(message = 
             paste0(v.name, " must not have missing values.\n"))
    }
  } 
  if(length.one) {
    if(length(v) != 1) {
      stop(message = paste0(v.name, " must be of length 1.\n"))
    }    
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
  
  if(!is.numeric(v)) {
    stop(message = 
           paste0(v.name, " must be numeric.\n"))
  }
  if(nonnegative) {
    if(any(v < 0)) {
      stop(message = 
             paste0(v.name, " must have nonnegative elements.\n"))
    }
  }
  if(positive) {
    if(any(v <= 0)) {
      stop(message = 
             paste0(v.name, " must have positive elements.\n"))
    }
  }
}

CheckUnitInterval <- function(v, v.name = NULL, 
  no.missing = T, length.one = T, equal = F) {
  # function that checks that a number is in (0,1) or [0,1].
  #  
  # Args: 
  #   v: an object to be tested if numeric.
  #   v.name: the name of v passed along from another function. 
  #   no.missing: an indicator that no missing values are allowed.
  #   length.one: an indicator that the argument is of length 1.
  #   equal: an indicator that the valid interval is [0, 1], not (0, 1).
  # Returns nothing.
  #
  # Error handling:
  v.name <- ifelse(is.null(v.name), deparse(substitute(v)), v.name)
  CheckCharacter(v.name)
  CheckBinary(no.missing)
  CheckBinary(length.one)
  CheckBinary(equal)
  
  if(sum(is.na(v)) == length(v)) {
    stop(message = paste0(v.name, " cannot contain only NA.\n"))
  }
  if(no.missing) {
    if(any(is.na(v))) {
      stop(message = 
             paste0(v.name, " must not have missing values.\n"))
    }
  }
  if(length.one) {
    if(length(v) != 1) {
      stop(message = paste0(v.name, " must be of length 1.\n"))
    }    
  }
  if(sum(is.na(v)) > 0) {
    v <- v[-which(is.na(v))]
  }
  if(!equal) {
    if(any(v <= 0)) {
      stop(message = paste0(v.name, " cannot be less than or equal to 0.\n"))
    }
    if(any(v >= 1)) {
      stop(message = paste0(v.name, " cannot be greater than or equal to 1.\n"))
    }    
  } else {
    if(any(v < 0)) {
      stop(message = paste0(v.name, " cannot be less than 0.\n"))
    }
    if(any(v > 1)) {
      stop(message = paste0(v.name, " cannot be greater than 1.\n"))
    }    
  }
}

CheckLocs <- function(locs) {
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
  CheckNumeric(locs, length.one = F)
  if(length(locs) != length(unique(locs))) {
    stop(message = "locs must be unique.\n")
  }
}

CheckInteriorKnots <- function(interior.knots) {
  # check that interior knot sequence is of correct form.
  # Args: 
  #   interior.knots: vector of knots that are in (0, 1). 
  # Returns nothing 
  # Error handling
  CheckUnitInterval(interior.knots, length.one = F)
  n.knots <- length(interior.knots)
  if(sum(order(interior.knots) == (1:n.knots)) < n.knots) {
    stop(message = "interior.knots must be ordered sequentially.\n")
  }
}

CheckCubicSplineKnots <- function(spline.knots) {
  # check that cubic spline knots is of correct form.
  # Args: 
  #   spline.knots: vector of knots. 
  # Returns nothing 
  # Error handling
  CheckUnitInterval(spline.knots, length.one = F, equal = T)
  n.knots <- length(spline.knots)
  if(n.knots < 9) {
    stop(message = "spline.knots must be a sequence of length at least 9.\n")
  }
  if(any(spline.knots[1:3] != 0)) {
    stop(message = "spline.knots must have first 3 elements equal to 0.\n")
  }
  if(any(spline.knots[(n.knots - 4):n.knots] != 1)) {
    stop(message = "spline.knots must have last 5 elements equal to 1.\n")
  }
  interior.knots <- spline.knots[4:(n.knots - 5)]
  CheckUnitInterval(interior.knots, length.one = F)
  if(any(order(interior.knots) != 1:length(interior.knots))) {
    stop(message = "spline.knots must be sorted in increasing order.\n")
  }
}

CheckInvertibleSymmetricMatrix <- function(v, v.name = NULL) {
  # Checks that a matrix is invertible and symmetric.
  #
  # Args: 
  #   v: matrix.
  # Returns nothing.
  # Error handling:
  v.name <- ifelse(is.null(v.name), deparse(substitute(v)), v.name)
  CheckCharacter(v.name)
  if(!is.matrix(v)) {
    stop(message = 
           paste0(v.name, " must be a matrix.\n"))
  }
  if(!isSymmetric(v)) {
    stop(message = 
           paste0(v.name, " must be symmetric.\n"))
  }
  flag <- 0
  v.inv <- 0
  suppressWarnings(
    tryCatch(
      v.inv <- chol2inv(
        chol(v)
      ), 
      error = function(e) {return(v.inv)}
    )
  )
  if (length(v.inv) == 1) {
    if (v.inv == 0) {
      flag <- 1
    }
  }
  if(flag == 1) {
    stop(message = 
           paste0(v.name, " is symmetric but not invertible.\n"))
  }
}

SetEqual <- function(x, c.indicator = F) {
  # creates a copy of numeric array x.
  # Args: 
  #    x: vector of reals.
  # Returns a copy of x.
  # Error handling:
  CheckNumeric(x, length.one = F)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    return(x)
  } else {
    length.x <- length(x)
    out <- .C("SetEqual",
              length_x = as.integer(length.x), 
              x = as.double(x), 
              copy_x  = as.double(rep(0, length.x))
    )
    x <- out$x
  }
  return(x)
}

Plot <- function(x, y) {
  # Plots vector y vs vector x.
  # Args: 
  #    x: vector of reals.
  #    y: vector of reals.
  # Returns nothing.
  # Error handling:
  CheckNumeric(x, length.one = F)
  CheckNumeric(y, length.one = F)
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

Mean <- function(x, c.indicator = F) {
  # Returns the mean of vector x.
  # Args: 
  #    x: vector of reals.
  #    c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   mean.x: the mean of x.
  # Error handling:
  CheckNumeric(x, length.one = F)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    mean.x <- mean(x)
  } else {
    out <- .C("Mean",
              length_x = as.integer(length(x)), 
              x = as.double(x), 
              mean_x  = as.double(0)
    )
    mean.x <- out$mean_x
  }
  return(mean.x)
}

Variance <- function(x, c.indicator = F) {
  # Returns the variance of vector x.
  # Args: 
  #    x: vector of reals.
  #    c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   mean.x: the mean of x.
  # Error handling:
  CheckNumeric(x, length.one = F)
  CheckBinary(c.indicator)
  length.x <- length(x)
  if(length.x == 1) {
    stop(message = "x must be of length greater than 1.\n")
  }
  if(!c.indicator) {
    var.x <- var(x)
  } else {
    out <- .C("Variance",
              length_x = as.integer(length(x)), 
              x = as.double(x), 
              var_x  = as.double(0)
    )
    var.x <- out$var_x
  }
  return(var.x)
}

FindIntegerMax <- function(x, c.indicator = F) {
  # finds the maximum value in a vector of integers.
  # Args:
  #   x: the vector of integers to search through.
  # Returns:
  #   max.x: the maximum value of x.
  # Error handling:
  CheckInteger(x)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    max.x <- max(x)
  } else {
    length.x <- length(x)
    out <- .C("FindIntegerMax",
              length_x = as.integer(length.x), 
              x = as.integer(x), 
              max_x = as.integer(0)
    )
    max.x <- out$max_x
  }
  return(max.x)
}

FindInterval <- function(tau, x, c.indicator = F) {
  # Finds which interval of x scalar tau is in
  #
  # Args: 
  #   tau: real scalar
  #   x: sorted vector of reals
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   bin: the bin (indexed from 0 to length(x) - 1) such that x[bin] <= tau < x[bin + 1].
  # For example: if x = 1:10, and tau = 1.5, then bin = 0;
  #              if x = 1:10, and tau = 9.5, then bin = 8;
  #              if x = 1:10, and tau = 100, then bin = 9;
  #  In R add 1 to the bins above.
  # Error handling:
  CheckNumeric(tau)
  CheckNumeric(x, length.one = F)
  length.x <- length(x)
  if(length.x != length(unique(x))) {
    stop(message = "values of x must be unique.\n")
  }
  if(any(order(x) != 1:length.x)) {
    stop(message = "x must be in ascending order.\n")
  }
  
  if(tau < x[1]) {
    stop(message = "tau must be greater than or equal to lowest value of x.\n")
  }
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    bin <- findInterval(tau, x)
  } else {
    out <- .C("FindInterval",
              length_x = as.integer(length.x), 
              x = as.double(x), 
              tau  = as.double(tau),
              bin = as.integer(1)
    )
    bin <- out$bin
  }
  return(bin)
}

MatrixMultiply <- function(a, b, trans.a = 0, trans.b = 0, c.indicator = F) {
  # Calculates the product of matrices a and b.
  #
  # Args: 
  #   a: matrix of reals.
  #   b: matrix of reals.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   c: product of a and b.
  #
  # Error handling:
  CheckNumeric(a, length.one = F)
  CheckNumeric(b, length.one = F)
  if(!is.matrix(a)) {
    stop(message = "a must be a matrix.\n")
  }
  if(!is.matrix(b)) {
    stop(message = "b must be a matrix.\n")
  }
  nrow.a <- nrow(a)
  ncol.a <- ncol(a)
  nrow.b <- nrow(b)
  ncol.b <- ncol(b)
  
  CheckBinary(trans.a)
  CheckBinary(trans.b)
  
  if(trans.a == 0 && trans.b == 0) {
    if(ncol.a != nrow.b) {
      stop(message = "number of columns of a must equal number of rows of b.\n")
    }
  }
  if(trans.a == 1 && trans.b == 0) {
    if(nrow.a != nrow.b) {
      stop(message = "number of rows of a must equal number of rows of b.\n")
    }
  }
  if(trans.a == 0 && trans.b == 1) {
    if(ncol.a != ncol.b) {
      stop(message = "number of columns of a must equal number of cols of b.\n")
    }
  }
  if(trans.a == 1 && trans.b == 1) {
    if(nrow.a != ncol.b) {
      stop(message = "number of rows of a must equal number of cols of b.\n")
    }
  }
  
  CheckBinary(c.indicator)
  if(!c.indicator) {
    if(trans.a) {
      a <- t(a)
    }
    if(trans.b) {
      b <- t(b)
    }
    c <- a %*% b
  } else {
    nrow.c <- ifelse(trans.a, ncol.a, nrow.a)
    ncol.c <- ifelse(trans.b, nrow.b, ncol.b)
    mid.c <- ifelse(trans.a, nrow.a, ncol.a)
    
    dim.c <- c(nrow.c, ncol.c)
    out <- .C("MatrixMultiply",
              t_a = as.integer(trans.a), 
              t_b = as.integer(trans.b), 
              nrow_a  = as.integer(nrow.a), 
              ncol_a  = as.integer(ncol.a), 
              nrow_b  = as.integer(nrow.b), 
              ncol_b  = as.integer(ncol.b), 
              a       = as.double(a),
              b       = as.double(b),
              c       = as.double(rep(0, nrow.c * ncol.c))
    )
    c <- array(out$c, dim = dim.c)
  }    
  return(c)
}

Kronecker <- function(a, b, c.indicator = F) {
  # Calculates the kronecker product of a and b.
  #
  # Args: 
  #   a: matrix of reals.
  #   b: matrix of reals.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   c: kronecker product of a and b.
  #
  # Error handling:
  CheckNumeric(a, length.one = F)
  CheckNumeric(b, length.one = F)
  if(!is.matrix(a)) {
    stop(message = "a must be a matrix.\n")
  }
  if(!is.matrix(b)) {
    stop(message = "b must be a matrix.\n")
  }
  CheckBinary(c.indicator)
  if(!c.indicator) {
    c <- kronecker(a, b)
  } else {
    nrow.a <- nrow(a)
    ncol.a <- ncol(a)
    nrow.b <- nrow(b)
    ncol.b <- ncol(b)
    
    out <- .C("Kronecker",
              nrow_a  = as.integer(nrow.a), 
              ncol_a  = as.integer(ncol.a), 
              nrow_b  = as.integer(nrow.b), 
              ncol_b  = as.integer(ncol.b), 
              a       = as.double(a),
              b       = as.double(b),
              c       = as.double(rep(0, nrow.a * nrow.b * ncol.a * ncol.b))
    )    
    c <- array(out$c, dim = c(nrow.a * nrow.b, ncol.a * ncol.b))
  }
  return(c)
}

QuadraticForm <- function(x, a, c.indicator = F) {
  # Finds qf = x'a x, where a is symmetric.
  #
  # Args: 
  #   x: vector.
  #   a: symmetric matrix that is length(x) x length(x).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   qf: the quadratic form.
  # Error handling:
  CheckNumeric(x, length.one = F)
  length.x <- length(x)
  CheckNumeric(a, length.one = F)
  if(!is.matrix(a)) {
    stop(message = "a must be a matrix.\n")
  }
  if(!isSymmetric(a)) {
    stop(message = "a must be symmetric.\n")
  }
  CheckBinary(c.indicator)
  if(!c.indicator) {
    qf <- t(x) %*% a %*% x
  } else {
    out <- .C("QuadraticForm",
              length_x = as.integer(length.x), 
              x = as.double(x), 
              a  = as.double(a),
              qf = as.double(0)
    )
    qf <- out$qf
  }
  return(qf)
}

QuadraticFormKronecker <- function(d.mat, e.mat, x, c.indicator = F) {
  # This function returns t(x) %*% kronecker(d.mat, e.mat) %*% x
  #   without calculating the kronecker product
  #
  # Args:
  #   d.mat : square matrix.
  #   e.mat : square matrix.
  #   x: vector of length nrow(d.mat) * nrow(e.mat).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns the quadratic form.
  # Error handling:
  CheckNumeric(d.mat, length.one = F)  
  if(!is.matrix(d.mat)) {
    stop(message = "d.mat must be a matrix.\n")
  }
  CheckNumeric(e.mat, length.one = F)  
  if(!is.matrix(e.mat)) {
    stop(message = "e.mat must be a matrix.\n")
  }
  nrow.d <- nrow(d.mat)
  if(nrow.d != ncol(d.mat)) {
    stop(message = "d.mat must be a square matrix.\n")
  }
  nrow.e <- nrow(e.mat)
  if(nrow.e != ncol(e.mat)) {
    stop(message = "e.mat must be a square matrix.\n")
  }
  CheckNumeric(x, length.one = F)
  if(length(x) != (nrow.d * nrow.e)) {
    stop(message = "x must be of length nrow(d) * nrow(e).\n")
  }
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    new.x <- matrix(x, nrow.e, nrow.d)
    matrix.product <- e.mat %*% new.x %*% d.mat
    qf <- sum(new.x * matrix.product)
    
#     qf <- 0
#     for(d in 1:nrow.d) {
#       qf <- qf + new.x[, d] %*% matrix.product[, d]
#     } 
  } else {
    out <- .C("QuadraticFormKronecker",
              nrow_d  = as.integer(nrow.d), 
              nrow_e  = as.integer(nrow.e), 
              d       = as.double(d.mat),
              e       = as.double(e.mat),
              x       = as.double(x),
              qf      = as.double(0)

    )
    qf <- out$qf
  }
  return(qf)  
}

CholMatrix <- function(v, c.indicator = F) {
  # Calculates the cholesky decompostion of v and the log determinant of v
  #
  # Args: 
  #   v: a positive definite symmetric matrix
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns: 
  #   v.chol: the lower triangular Cholesky decomposition of v.
  #   log.det: the log determinant of v.
  # Some error handling done in C.
  CheckNumeric(v, length.one = F)
  if(!is.matrix(v)) {
    stop(message = "v must be a matrix.\n")
  }
  
  CheckBinary(c.indicator)

  if(!c.indicator) {
    v.chol <- t(chol(v))
    log.det <- 2 * sum(log(diag(v.chol)))
  }

  else {
    nrow.v <- nrow(v)
    out <- .C("CholMatrix",
              nrow_v   = as.integer(nrow.v),
              v        = as.double(v),            
              v_chol   = as.double(rep(0, nrow.v * nrow.v)),
              log_det   = as.double(0)
    )
    v.chol <- array(out$v_chol, dim = c(nrow.v, nrow.v))
    log.det <- out$log_det
  }

  return(
    list(v.chol = v.chol, log.det = log.det)
  )
}

InvertSymmetricMatrix <- function(v, c.indicator = F) {
  # Calculates the inverse of v and the log determinant of the inverse.
  #
  # Args: 
  #   v: an invertible symmetric matrix.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns: 
  #   v.inv: the inverse of v.
  #   log.det.v.inv: the log determinant of the inverse.
  # Some error handling done in C.
  
  CheckNumeric(v, length.one = F)
  if(!is.matrix(v)) {
    stop(message = "v must be a matrix.\n")
  }
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    chol.v <- chol(v) #upper triangular form
    log.det <- -2 * sum(log(diag(chol.v)))
    v.inv <- chol2inv(chol.v)
  } else {
    nrow.v <- nrow(v)
    out <- .C("InvertSymmetricMatrix",
              nrow_v   = as.integer(nrow.v),
              v        = as.double(v),            
              v_inv   = as.double(rep(0, nrow.v * nrow.v)),
              log_det   = as.double(0)
    )
    v.inv <- array(out$v_inv, dim = c(nrow.v, nrow.v)) 
    log.det <- out$log_det    
  }
  list(v.inv = v.inv, log.det = log.det)
}

MakeAbsoluteDistance <- function(timepoints, c.indicator = F) {
  # Constructs absolute distance matrix using unique vector of timepoints. 
  #
  # Args: 
  #   timepoints: vector of unique timepoints. 
  # Returns: 
  #   distance.matrix: autoregressive distance matrix.  
  #   Error handling:  
  if(FALSE %in% (timepoints == sort(timepoints))) {
    stop(message = "timepoints must be ordered sequentially.\n")
  }
  if(length(timepoints) != length(unique(timepoints))) {
    stop(message = "timepoints must contain unique values.\n")
  }
  if(!c.indicator) {
    distance.matrix <- abs(outer(timepoints, timepoints, "-"))  
  } else {
    n.timepoints <- length(timepoints)
    out <- .C("MakeAbsoluteDistance",
              n_timepoints = as.integer(n.timepoints),
              timepoints   = as.double(timepoints),
              distance_matrix = as.double(rep(0, n.timepoints * n.timepoints))
    )
    distance.matrix <- matrix(out$distance_matrix, n.timepoints, n.timepoints)
  }
  return(distance.matrix)
}

MakeAutoregressivePrecision <- function(rho, distance.matrix, c.indicator = F){
  # Constructs precision matrix of AR-1 process. 
  #
  # Args: 
  #   rho: scalar in (0, 1).    
  #   distance.matrix: matrix indicating the distance between timepoints. 
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns: 
  #   prec.ar: autoregressive precision matrix.  
  #   log.det: log determinant of the precision matrix.
  # Error handling:
  if(length(rho) > 1){
    stop(message = "rho must be of length 1.\n")
  }
  CheckUnitInterval(rho)
  if(!is.matrix(distance.matrix)) {
    stop(message = "distance.matrix must be a matrix.\n")
  }
  CheckNumeric(distance.matrix, length.one = F)
  if(!isSymmetric(distance.matrix)) {
    stop(message = "distance.matrix must be symmetric.\n")
  } 
  if(any(diag(distance.matrix) != 0)) {
    stop(message = "distance.matrix must have a diagonal of 0.\n")
  }
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    ar <- rho ^ distance.matrix
    out <- InvertSymmetricMatrix(ar)
    prec.ar <- out$v.inv
    log.det <- out$log.det
  } else {
    n.timepoints <- nrow(distance.matrix)
    out <- .C("MakeAutoregressivePrecision",
              n_timepoints  = as.integer(n.timepoints), 
              rho           = as.double(rho),
              dist          = as.double(distance.matrix),
              prec_ar       = as.double(rep(0, n.timepoints * n.timepoints)),
              log_det       = as.double(0)
    )
    prec.ar <- array(out$prec_ar, dim = c(n.timepoints, n.timepoints))
    log.det <- out$log_det
  }
  return(list(prec.ar = prec.ar, log.det = log.det))
}

MakeSpatialDistance <- function(locs) {
  # Constructs Euclidean distance matrix using unique vector of locations. 
  #
  # Args: 
  #   locs: matrix with 2 columns of unique locations. 
  # Returns: 
  #   spatial.dist: Euclidean distance matrix.  
  # Error handling:  
  CheckLocs(locs)
  
  n.locs <- nrow(locs)
  lon.dist <- outer(locs[,1], locs[,1], "-")
  lat.dist <- outer(locs[,2], locs[,2], "-")
  spatial.dist <- sqrt(lon.dist ^ 2 + lat.dist ^ 2)
  return(spatial.dist)
}

MakeSpatialPrecision <-function(rho, spatial.dist, c.indicator = F){
  # Constructs exponential spatial precision matrix. 
  #
  # Args: 
  #   rho: range parameter (positive scalar).  
  #   spatial.dist: euclidean distance matrix.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns: 
  #   prec.sp: spatial precision matrix.  
  #   log.det: log determinant of the precision matrix.
  # 
  # Error handling:
  if(length(rho) > 1){
    stop(message = "rho must be of length 1.\n")
  }
  CheckNumeric(rho, length.one = T, positive = T)
  if(!is.matrix(spatial.dist)) {
    stop(message = "spatial.dist must be a matrix.\n")
  }
  CheckNumeric(spatial.dist, length.one = F)
  if(!isSymmetric(spatial.dist)) {
    stop(message = "spatial.dist must be symmetric.\n")
  } 
  if(any(diag(spatial.dist) != 0)) {
    stop(message = "spatial.dist must have a diagonal of 0.\n")
  }
  CheckBinary(c.indicator)
  if(!c.indicator) {
    sp <- exp(-spatial.dist / rho)
    out <- InvertSymmetricMatrix(sp)
    prec.sp <- out$v.inv
    log.det <- out$log.det
  } else {
    n.locs <- nrow(spatial.dist)
    out <- .C("MakeSpatialPrecision",
              n_locs    = as.integer(n.locs), 
              rho       = as.double(rho),
              dist      = as.double(spatial.dist),
              prec_sp   = as.double(rep(0, n.locs * n.locs)),
              log_det   = as.double(0)
    )
    prec.sp <- array(out$prec_sp, dim = c(n.locs, n.locs))
    log.det <- out$log_det
  }
  return(list(prec.sp = prec.sp, log.det = log.det))
}

MakeSplineKnots <- function(interior.knots, spline.df = 3, c.indicator = F) {
  #constructs knot sequence for spline basis.
  # Args: 
  #   interior.knots = vector of knots that are in (0, 1). 
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns: 
  #   spline.knots: vector of knots for the spline basis.
  # Error handling:
  CheckInteriorKnots(interior.knots)
  CheckInteger(spline.df, positive = T)
  CheckBinary(c.indicator)

  if(!c.indicator) {
    spline.knots <- c(rep(0, spline.df), interior.knots, rep(1, spline.df + 2))    
  } else {
    interior.knots.length <- length(interior.knots)
    out <- .C("MakeSplineKnots",
              interior_knots_length = as.integer(interior.knots.length), 
              spline_df             = as.integer(spline.df), 
              interior_knots        = as.double(interior.knots), 
              spline_knots          = as.double(rep(0, 2 * spline.df + interior.knots.length + 2))
    )
    spline.knots <- out$spline_knots
  }    
  return(spline.knots)
}

MSpline <- function(tau, spline.df, spline.knots, m, c.indicator = F) {
  # Given knots spline.knots, evaluates the mth m-spline of  
  #   spline.df degrees of freedom at tau.
  # Args:
  #   tau: scalar in (0, 1).
  #   spline.df: degrees of freedom for splines.
  #   spline.knots: knot sequence.
  #   m: index for the spline.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interior.knots and spline.df degrees of freedom, 
  #       the mth m-spline evaluated at tau.
  # Error handling
  CheckUnitInterval(tau)
  CheckInteger(spline.df, positive = T)
  CheckNumeric(spline.knots, length.one = F, nonnegative = T)
  if(min(spline.knots) != 0) {
    stop(message = "spline.knots first element must be 0.\n")
  }
  if(max(spline.knots) != 1) {
    stop(message = "spline.knots greatest element must be 1.\n")
  }
  if(any(order(spline.knots) != (1:length(spline.knots)))) {
    stop(message = "spline.knots must be in increasing order.\n")
  }
  
  CheckInteger(m, positive = T)
  if(m  > length(spline.knots) - spline.df + 1) {
    stop(message = "m cannot be greater than length(spline.knots) - spline.df + 1.\n")
  }

  CheckBinary(c.indicator)

  if(!c.indicator) {
    if(tau < spline.knots[m] || tau >= spline.knots[m + spline.df]) {
      v <- 0
    } else if (spline.df == 1) {
      v <- 1 / (spline.knots[m + 1] - spline.knots[m])
    } else {
      d1 <- tau - spline.knots[m]
      d2 <- spline.knots[m + spline.df] - tau
      v <- d1 * MSpline(tau, spline.df - 1, spline.knots, m) + 
           d2 * MSpline(tau, spline.df - 1, spline.knots, m + 1)
      v <- v * spline.df 
      v <- v / ((spline.df - 1) * (spline.knots[m + spline.df] - spline.knots[m]))
    }
  } else {
    out <- .C("MSpline",
              tau                 = as.double(tau),            
              spline_df           = as.integer(spline.df),            
              spline_knots        = as.double(spline.knots),            
              m                   = as.integer(m - 1),
              v                   = as.double(0)            
    )    
    v <- out$v
  }
  return(v)
}

ISpline <- function(tau, spline.df, spline.knots, m, c.indicator = F) {
  # Given knots spline.knots, evaluates the mth i-spline of  
  #   spline.df degrees of freedom at tau.
  # Args:
  #   tau: scalar in (0, 1).
  #   spline.df: degrees of freedom for splines.
  #   spline.knots: knot sequence.
  #   m: index for the spline.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interior.knots and spline.df degrees of freedom, 
  #       the mth i-spline evaluated at tau.
  # Error handling
  CheckUnitInterval(tau)
  CheckInteger(spline.df, positive = T)
  
  CheckNumeric(spline.knots, length.one = F, nonnegative = T)
  if(min(spline.knots) != 0) {
    stop(message = "spline.knots first element must be 0.\n")
  }
  if(max(spline.knots) != 1) {
    stop(message = "spline.knots greatest element must be 1.\n")
  }
  if(any(order(spline.knots) != (1:length(spline.knots)))) {
    stop(message = "spline.knots must be in increasing order.\n")
  }
  
  CheckInteger(m, positive = T)
  if(m  > length(spline.knots) - spline.df) {
    stop(message = "m cannot be greater than length(spline.knots) - 3.\n")
  }
  
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    bin <- findInterval(tau, spline.knots)  
      if (bin < m){ 
        v<- 0
      } else if ((bin - spline.df + 1) > m ){
        v <- 1
      } else {
        v <- 0
        for (l in m : (bin + 1)){
          v <- v + MSpline(tau, spline.df + 1, spline.knots, l) * 
            (spline.knots[l + spline.df + 1] - spline.knots[l]) / 
            (spline.df + 1)
        }
      }
    } else {
      spline.knots.length <- length(spline.knots)
      out <- .C("ISpline",
                tau                 = as.double(tau),            
                spline_df           = as.integer(spline.df),            
                spline_knots_length = as.integer(spline.knots.length),            
                spline_knots        = as.double(spline.knots),            
                m                   = as.integer(m - 1),
                v                   = as.double(0)            
      )    
      v <- out$v
  }
  return(v)
}

CubicMSpline <- function(tau, interior.knots, m, c.indicator = F) {
  # Given knots spline.knots, evaluates the mth m-spline of  
  #   3 degrees of freedom at tau.
  # Args:
  #   tau: scalar in (0, 1).
  #   spline.df: degrees of freedom for splines (should be 3).
  #   interior.knots: interior knot sequence.
  #   m: index for the spline.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interior.knots and spline.df degrees of freedom, 
  #       the mth m-spline evaluated at tau
  # Error handling
  CheckUnitInterval(tau)
  CheckInteriorKnots(interior.knots)
  CheckInteger(m, positive = T)
  spline.knots <- MakeSplineKnots(interior.knots, spline.df = 3, c.indicator = F)
  
  if(m  > length(spline.knots) - 3) {
    stop(message = "m cannot be greater than length(spline.knots) - 3.\n")
  }
  CheckBinary(c.indicator)
  spline.df <- 3
  if(!c.indicator) {
    bin <- findInterval(tau, spline.knots)
    if (bin < m) {
      v <- 0
    } else if (
      (bin - spline.df + 1) > m) {
      v <- 0
    } else if (bin == m) {
      v <- 3 * (tau - spline.knots[m]) ^ 2 /
        ((spline.knots[m + 1] - spline.knots[m]) * 
           (spline.knots[m + 2] - spline.knots[m]) * 
           (spline.knots[m + 3] - spline.knots[m])
        )
    } else if (bin == (m + 1)) {
      i1 <-  1 / (spline.knots[m + 2] - spline.knots[m + 1])
      i2 <-  (2 * (spline.knots[m + 3] - spline.knots[m + 1])) / 
        ((spline.knots[m + 3] -spline.knots[m + 1]) * 
           (spline.knots[m + 2] - spline.knots[m + 1])
        )
      i3 <-  -3 *(
        (spline.knots[m + 3] - tau) ^ 2/
          ((spline.knots[m + 3] - spline.knots[m + 1]) * 
             (spline.knots[m + 3] - spline.knots[m]) * 
             (spline.knots[m + 2] - spline.knots[m + 1])) + 
          (tau - spline.knots[m]) ^ 2/
          (
            (spline.knots[m + 3] - spline.knots[m]) * 
              (spline.knots[m + 2] - spline.knots[m]) * 
              (spline.knots[m + 2] - spline.knots[m + 1])
          )      
      )
      v <- i1 + i2 + i3
    } else{
      v <- 3 * (spline.knots[m + 3] - tau) ^ 2 / 
        ((spline.knots[m + 3] - spline.knots[m + 2]) * 
           (spline.knots[m + 3] - spline.knots[m + 1]) * 
           (spline.knots[m + 3] - spline.knots[m])
        )
    }  
  } else {
    out <- .C("CubicMSpline",
              tau                 = as.double(tau),            
              spline_knots_length = as.integer(length(spline.knots)),            
              spline_knots        = as.double(spline.knots),            
              m                   = as.integer(m - 1),
              bin                 = as.integer(0),
              v                   = as.double(0)            
    )    
    v <- out$v
  }
  return(v)
}

CubicISpline <- function(tau, interior.knots, m, c.indicator = F){
  # Given knots spline.knots, evaluates the mth i-spline 
  #  of 3 degrees of freedom at tau.
  #   
  # Args:
  #   tau: scalar in (0, 1).
  #   interior.knots: interior knot sequence.
  #   m: index for the spline.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interior.knots and spline.df degrees of freedom, 
  #       the mth i-spline evaluated at tau.
  # Error handling
  CheckUnitInterval(tau)
  CheckInteriorKnots(interior.knots)
  CheckInteger(m, positive = T)
  spline.knots <- MakeSplineKnots(interior.knots, spline.df = 3, c.indicator = F)
  if(m  > length(spline.knots) - 3) {
    stop(message = "m cannot be greater than length(spline.knots) - 3.\n")
  }
  CheckBinary(c.indicator)
  spline.df <- 3
  if(!c.indicator) {
    bin <- findInterval(tau,spline.knots)
    if(bin < m){
      v <- 0
    }
    else if ((bin - spline.df + 1) > m) {
      v <- 1
    }
    else if (bin == m){
      v <- ((tau - spline.knots[m]) ^ 3) /
        (
          (spline.knots[m+1]  - spline.knots[m]) * 
            (spline.knots[m+2]  - spline.knots[m]) * 
            (spline.knots[m+3]  - spline.knots[m])
        )
    }
    else if(bin == (m + 1)){
      i1 <- (tau - spline.knots[m]) / (spline.knots[m + 2] - spline.knots[m + 1])
      i2 <- ((tau - spline.knots[m + 1]) ^ 2 - 
               (spline.knots[m + 3]-tau) ^ 2) /
        ((spline.knots[m + 3] - spline.knots[m + 1]) * (spline.knots[m + 2] - spline.knots[m + 1]))
      i3 <- (spline.knots[m + 3] - tau) ^ 3 / 
        ((spline.knots[m + 3] - spline.knots[m + 1]) * 
           (spline.knots[m + 3] - spline.knots[m]) * 
           (spline.knots[m + 2]-spline.knots[m + 1])) - 
        (tau - spline.knots[m]) ^ 3 / 
        ((spline.knots[m + 3] - spline.knots[m]) * 
           (spline.knots[m + 2] - spline.knots[m]) * 
           (spline.knots[m+2]-spline.knots[m+1])
        )
      v <- i1 + i2 + i3
    }
    else{
      v <- 1 - (spline.knots[m + 3] - tau) ^ 3 / 
        ((spline.knots[m + 3] - spline.knots[m + 2]) *
           (spline.knots[m + 3] - spline.knots[m + 1]) *
           (spline.knots[m + 3] - spline.knots[m]))
    }
  } else {
    out <- .C("CubicISpline",
              tau                 = as.double(tau),            
              spline_knots_length = as.integer(length(spline.knots)),            
              spline_knots        = as.double(spline.knots),            
              m                   = as.integer(m - 1),
              bin                 = as.integer(0),
              v                   = as.double(0)            
    )    
    v <- out$v
  }
  return(v)
}

MakeCubicSplineCoefficients <- function(interior.knots, c.indicator = F) {
  # Constructs array containing coefficients for cubic spline calculation.
  # Args: 
  #   interior.knots: a sequence of knots in (0, 1).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   spline.coefs: an array of dimension 4 x n.bf x (n.bf - 1)
  #     containing the coefficients for cubic spline calculation.  
  #     spline.coefs[1, , ] contains the cubic coefficients. 
  #     spline.coefs[2, , ] contains the quadratic coefficients.
  #     spline.coefs[3, , ] contains the linear coefficients.
  #     spline.coefs[4, , ] contains the constant coefficients.
  #     The second dimension indexes the basis function.
  #     The third dimension indexes the bin.  
  #     Given that tau is in the mth bin, 
  #       I_l(tau) =  spline.coefs[1,l, m] * tau ^ 3 
  #           + spline.coefs[2,l, m] * tau ^ 2
  #           + spline.coefs[3,l, m] * tau 
  #           + spline.coefs[4,l, m]. 
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
  # d.a, d.b, d.c, and d.d are (n.bf - 1) x (n.bf - 1) matrices.  
  # Row indexes basis function (spline, that is).  
  # column indexes bin.
  # d.a represents the coefficients for the cubic coefficient terms.
  # d.b represents the quadratic coefficient.
  # d.c represents the linear coefficient.
  # d.d represents the constant coefficient.
  #Error Handling
  CheckInteriorKnots(interior.knots)
  CheckBinary(c.indicator)
  spline.knots <- MakeSplineKnots(interior.knots)  
  n.bf <- length(spline.knots) - 4
  if(!c.indicator) {
    a <- matrix(0, nrow = length(spline.knots), 1)
    b <- matrix(0, nrow = length(spline.knots), 4)
    c <- matrix(0, nrow = length(spline.knots), 1)
    for(m in 3:(n.bf - 1)) {
      a[m] <- 1 / ((spline.knots[m + 1] - spline.knots[m]) * 
                     (spline.knots[m + 2] - spline.knots[m]) * 
                     (spline.knots[m + 3] - spline.knots[m]) 
      )
    }
    for(bin in 3:(n.bf - 1)) {
      m <- bin - 1
      b[bin,1] <- 1 / (spline.knots[m + 2] - spline.knots[m + 1])
      b[bin,2] <- 1 /((spline.knots[m + 3] - spline.knots[m + 1])  * (spline.knots[m + 2] - spline.knots[m + 1]))
      b[bin,3] <- - 1/((spline.knots[m + 3] - spline.knots[m + 1]) * (spline.knots[m + 3] - spline.knots[m]) * (spline.knots[m + 2] - spline.knots[m + 1]))
      b[bin,4] <- - 1/((spline.knots[m + 3] - spline.knots[m]) * (spline.knots[m + 2] - spline.knots[m]) * (spline.knots[m + 2] - spline.knots[m + 1]))
    }
    for(bin in 3:(n.bf - 1)) {
      m <- bin - 2
      c[bin]   <-  1 / ((spline.knots[m + 3] - spline.knots[m + 2]) * (spline.knots[m + 3] - spline.knots[m + 1]) * (spline.knots[m+3] - spline.knots[m]))
    }
    
    d.a <- d.b <- d.c <- d.d <- matrix(0, nrow = n.bf - 1, ncol = n.bf - 1)
    for(l in 1:n.bf) { #indexes basis function
      for(m in 3: (n.bf - 1)) { #indexes bin
        if (m > (l + 2)) {d.d[l,m] <- 1}
        if (l == m) {
          d.a[l, m] <- a[m] 
          d.b[l, m] <- - 3 * a[m] * spline.knots[m]
          d.c[l, m] <-   3 * a[m] * spline.knots[m] ^ 2  
          d.d[l, m] <- - a[m] * spline.knots[m] ^ 3
        }
        if ((l + 1) == m) {
          d.a[l, m] <- b[m, 3] + b[m, 4] 
          d.b[l, m] <- - 3 * b[m, 3] * spline.knots[m + 2] - 3 * b[m, 4] * spline.knots[m - 1] 
          d.c[l, m] <-   3 * b[m, 3] * (spline.knots[m + 2] ^ 2) + 3 * b[m, 4] * (spline.knots[m - 1] ^ 2) + b[m,1] + 2 * b[m,2] * (spline.knots[m + 2] - spline.knots[m])   
          d.d[l, m] <- - b[m, 3] * spline.knots[m + 2] ^ 3 - b[m, 4] * (spline.knots[m - 1] ^ 3)  - b[m, 1] * spline.knots[m - 1]  + b[m, 2] * (spline.knots[m] ^ 2 - spline.knots[m + 2] ^ 2 )
        }
        if ((l + 2) == m) {
          d.a[l, m] <- c[m] 
          d.b[l, m] <- -3 * c[m] * spline.knots[m + 1]
          d.c[l, m] <- 3 * c[m] * spline.knots[m + 1]^2 
          d.d[l, m] <- -c[m] * spline.knots[m + 1] ^ 3 + 1
        }
      }
    }
    
    spline.coefs <- array(0, dim = c(4, n.bf, n.bf - 1))
    spline.coefs[1, 2 : n.bf,] <- d.a
    spline.coefs[2, 2 : n.bf,] <- d.b
    spline.coefs[3, 2 : n.bf,] <- d.c
    spline.coefs[4, 2 : n.bf,] <- d.d
    
    spline.coefs[4, 1, ] <- 1
  } else {
    out <- .C("MakeCubicSplineCoefficients",
              n_bf          = as.integer(n.bf),            
              spline_knots  = as.double(spline.knots),
              spline_coefs  = as.double(rep(0, 4 * n.bf * (n.bf - 1)))
    )       
    spline.coefs <- array(out$spline_coefs, dim = c(4, n.bf, n.bf - 1))
  }
  return(spline.coefs)
}

CubicMSpline2 <- function(tau, spline.knots, spline.coefs, m,  c.indicator = F) {
  # Given knots spline.knots and spline.coefs, evaluates the mth m-spline 
  #  of 3 degrees of freedom at tau.
  #   
  # Args:
  #   tau: scalar in (0, 1).
  #   spline.knots: knot sequence.
  #   spline.coefs: coefficients used to calculate the cubic spline.
  #      Assumed to be correctly calculated inside of this function.
  #   m: index for the spline.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   v : Given interior.knots and spline.df degrees of freedom, 
  #       the mth m-spline evaluated at tau.
  # Error handling
  CheckUnitInterval(tau)
  CheckCubicSplineKnots(spline.knots)
  n.bf <- length(spline.knots) - 4
  CheckNumeric(spline.coefs, length.one = F)
  if(any(dim(spline.coefs) != c(4, n.bf, n.bf - 1))) {
    stop(message = "spline.coefs must be a 4 x n.bf x n.bf - 1 array.\n")
  }
  CheckInteger(m, positive = T)
  if(m  > length(spline.knots) - 3) {
    stop(message = "m cannot be greater than length(spline.knots) - 3.\n")
  }
  CheckBinary(c.indicator)
  spline.df <- 3
  if(!c.indicator) {
    bin <- findInterval(tau, spline.knots)
    a <- spline.coefs[1, m + 1, bin]
    b <- spline.coefs[2, m + 1, bin]
    c <- spline.coefs[3, m + 1, bin]
    d <- spline.coefs[4, m + 1, bin]
    v <- 3 * a * tau ^ 2 + 2 * b * tau  + c
  } else {
    out <- .C("CubicMSpline2",
              tau                 = as.double(tau),
              spline_knots_length = as.integer(length(spline.knots)),            
              spline_knots        = as.double(spline.knots),
              spline_coefs        = as.double(spline.coefs),
              m                   = as.integer(m - 1),
              bin                 = as.integer(2),
              v                   = as.double(0)
    )       
    v <- out$v
  }
  return(v)
}    

CubicISpline2 <- function(tau, spline.knots, spline.coefs, m,  c.indicator = F) {
# Given knots spline.knots and spline.coefs, evaluates the mth i-spline 
#  of 3 degrees of freedom at tau.
#   
# Args:
#   tau: scalar in (0, 1).
#   spline.knots: knot sequence.
#   spline.coefs: coefficients used to calculate the cubic spline.
#      Assumed to be correctly calculated inside of this function.
#   m: index for the spline.
#   c.indicator: an indicator if the computation should be performed in C.
# Returns:
#   v : Given interior.knots and spline.df degrees of freedom, 
#       the mth i-spline evaluated at tau.
# Error handling:
  CheckUnitInterval(tau)
  CheckCubicSplineKnots(spline.knots)
  n.bf <- length(spline.knots) - 4
  CheckNumeric(spline.coefs, length.one = F)
  if(any(dim(spline.coefs) != c(4, n.bf, n.bf - 1))) {
    stop(message = "spline.coefs must be a 4 x n.bf x n.bf - 1 array.\n")
  }
  CheckInteger(m, positive = T)
  if(m  > length(spline.knots) - 3) {
    stop(message = "m cannot be greater than length(spline.knots) - 3.\n")
  }
  CheckBinary(c.indicator)
  spline.df <- 3
  if(!c.indicator) {
    bin <- findInterval(tau, spline.knots)
    a <- spline.coefs[1, m + 1, bin]
    b <- spline.coefs[2, m + 1, bin]
    c <- spline.coefs[3, m + 1, bin]
    d <- spline.coefs[4, m + 1, bin]
    v <- a * tau ^ 3 + b * tau ^ 2 + c * tau + d
  } else {
    out <- .C("CubicISpline2",
              tau                 = as.double(tau),
              spline_knots_length = as.integer(length(spline.knots)),            
              spline_knots        = as.double(spline.knots),
              spline_coefs        = as.double(spline.coefs),
              m                   = as.integer(m - 1),
              bin                 = as.integer(2),
              v                   = as.double(0)
    )       
    v <- out$v
  }
  return(v)
}    

MakeStickBreakingWeights <- function(v, c.indicator = F){
  # Transforms a vector of beta random variables and a 1 at the end into  
  #  latent class weights.
  #   
  # Args:
  #   v: vector containing beta random variables, 
  #      except for the final element which is 1.
  #      
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   nu : Latent class weights. 
  # Error handling:
  CheckNumeric(v, length.one = F)
  length.v <- length(v)
  if(v[length.v] != 1) {
    stop(message = "The last element of v must be 1.\n")
  }
  if(length.v > 1) {
    CheckUnitInterval(v[1 : (length.v - 1)], length.one = F)  
  }
  CheckBinary(c.indicator)
  if(!c.indicator) {
    if(length.v == 1) {
      nu <- v
    } else {
      nu <- v
      nu[2 : length.v] <- nu[2 : length.v] * cumprod(1 - v[1:(length.v - 1)])
    }
  } else {
    out <- .C("MakeStickBreakingWeights",
              length_v  = as.integer(length.v),
              v         = as.double(v),
              nu        = as.double(rep(0, length.v))            
    )       
    nu <- out$nu
  }
  return(nu)
}

LogDNorm <- function(y, loc, scale, shape = 0, c.indicator = F) {
  # Finds fY(y), where Y is normally distributed with 
  ## location scale parameters loc/scale
  #
  # Args: 
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter (ignored).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    d.y <- dnorm(y, loc, sd = scale, log = T)
  } else {
    out <- .C("LogDNorm",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              d_y = as.double(0)
    )    
    d.y <- out$d_y
  }
  return(d.y)
}

LogDT <- function(y, loc, scale, shape, c.indicator = F) {
  # Finds fY(y), where Y is student's t distributed with 
  ## location/scale/df parameters loc/scale/shape
  #
  # Args: 
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    z <- (y - loc) / scale
    d.y <- dt(z, df = shape, log = T) - log(scale)
  } else {
    out <- .C("LogDT",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              d_y = as.double(0)
    )
    d.y <- out$d_y
  }
  return(d.y)
}

LogDLogistic <- function(y, loc, scale, shape = 0, c.indicator = F) {
  # Finds fY(y), where Y is logistically distributed with 
  ## location scale parameters loc/scale
  #
  # Args: 
  #   y: scalar
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter (ignored).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    d.y <- dlogis(y, loc, scale, log = T)
  } else {
    out <- .C("LogDLogistic",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              d_y = as.double(0)
    )
    d.y <- out$d_y    
  }
  return(d.y)
}

LogDAsymmetricLaplace <- function(y, loc, scale, shape, c.indicator = F) {
  # Finds fY(y), where Y is asymmetric Lapace distributed with 
  ## location/scale/shape parameters loc/scale/shape
  #
  # Args: 
  #   c.indicator: an indicator if the computation should be performed in C.
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter in (0, 1).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    d.y <- log(dalap(y, loc, scale = scale, tau = shape))
  } else {
    out <- .C("LogDAsymmetricLaplace",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              d_y = as.double(0)
    )
    d.y <- out$d_y
  }
  return(d.y)
}

LogDWeibull <- function(y, loc = 0, scale, shape, c.indicator = F) {
  # Finds fY(y), where Y is Weibull distributed with 
  ## scale/shape parameters scale/shape.
  #
  # Args: 
  #   y: scalar.
  #   loc: location parameter (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    d.y <- dweibull(y, shape, scale, log = T)
  } else {
    out <- .C("LogDWeibull",
              y = as.double(y), 
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              d_y = as.double(0)
    )
    d.y <- out$d_y
  }
  return(d.y)
}

LogDGamma <- function(y, loc = 0, scale, shape, c.indicator = F) {
  # Finds fY(y), where Y is gamma distributed with 
  ## scale/shape parameters scale/shape. 
  #
  # Args: 
  #   y: scalar.
  #   loc: positive scale parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    d.y <- dgamma(y, shape = shape, scale = scale, log = T)
  } else {
    out <- .C("LogDGamma",
              y = as.double(y), 
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              d_y = as.double(0)
    )
    d.y <- out$d_y
  }
  return(d.y)
}

LogDBeta <- function(y, shape1, shape2, c.indicator = F) {
  # Finds fY(y), where Y is beta distributed with 
  ## parameters shape1 and shape2. 
  #
  # Args: 
  #   y: scalar.
  #   shape1: positive first shape parameter.
  #   shape1: positive second shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   log(fY(y)).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(shape1, positive = T)
  CheckNumeric(shape2, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    d.y <- dbeta(y, shape1, shape2, log = T)
  } else {
    out <- .C("LogDBeta",
              y = as.double(y), 
              shape1  = as.double(shape1),
              shape2  = as.double(shape2),
              d_y = as.double(0)
    )
    d.y <- out$d_y
  }
  return(d.y)
}

LogDInverseWishart <- function(df, sai, x, c.indicator = F) {
  # evaluates the log density of the inverse wishart distribution evaluated at x
  #   with scale matrix sai and df degrees of freedom.
  # Args:
  #   df = degrees of freedom.
  #   sai = scale matrix.
  #   x = positive definite matrix argument.
  # Returns:
  #  d.y: log density.
  # Error handling:
  CheckNumeric(df, positive = T)
  CheckInvertibleSymmetricMatrix(sai)
  CheckInvertibleSymmetricMatrix(x)
  p <- nrow(x)
  if(nrow(sai) != p) {
    stop(message = "sai must have the same number of rows as x.\n")
  }
  if(df <= (p - 1)) {
    stop(message = "df must be greater than p - 1.\n")
  }
  CheckBinary(c.indicator)
  if(!c.indicator) {
    d.y <- log(diwish(x, df, sai))
  } else {
    out <- .C("LogDInverseWishart",
              p = as.integer(p), 
              df = as.double(df),
              sai = as.double(sai),
              x = as.double(x),
              d_y = as.double(0)
    )
    d.y <- out$d_y
  }
  return(d.y) 
}    

PNorm <- function(y, loc, scale, shape = 0, c.indicator = F) {
  # Finds P(Y <= y), where Y is normally distributed with 
  ## location scale parameters loc/scale.
  #
  # Args: 
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    p.y <- pnorm(y, loc, sd = scale)
  } else {
    out <- .C("PNorm",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              p_y = as.double(0)
    )
    p.y <- out$p_y
  }
  return(p.y)
}

PT <- function(y, loc, scale, shape, c.indicator = F) {
  # Finds P(Y <= y), where Y is student's t distributed with 
  ## location/scale/df parameters loc/scale/shape.
  #
  # Args: 
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    p.y <- pt((y - loc) / scale, df = shape)
  } else {
    out <- .C("PT",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              p_y = as.double(0)
    )
    p.y <- out$p_y
  }
  return(p.y)
}

PLogistic <- function(y, loc, scale, shape = 0, c.indicator = F) {
  # Finds P(Y <= y), where Y is logistically distributed with 
  ## location scale parameters loc/scale.
  #
  # Args: 
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    p.y <- plogis(y, loc, scale)
  } else {
    out <- .C("PLogistic",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              p_y = as.double(0)
    )
    p.y <- out$p_y
  }
  return(p.y)
}

PAsymmetricLaplace <- function(y, loc, scale, shape, c.indicator = F) {
  # Finds P(Y <= y), where Y is asymmetric Lapace distributed with 
  ## location/scale/shape parameters loc/scale/shape.
  #
  # Args: 
  #   y: scalar.
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter in (0, 1).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    p.y <- palap(y, loc, scale, tau = shape)
  } else {
    out <- .C("PAsymmetricLaplace",
              y = as.double(y), 
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              p_y = as.double(0)
    )
    p.y <- out$p_y
  }
  return(p.y)
}

PWeibull <- function(y, loc = 0, scale, shape, c.indicator = F) {
  # Finds P(Y <= y), where Y is Weibull distributed with 
  ## scale/shape parameters scale/shape,
  #
  # Args: 
  #   y: scalar.
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y)
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    p.y <- pweibull(y, shape, scale)
  } else {
    out <- .C("PWeibull",
              y = as.double(y), 
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              p_y = as.double(0)
    )
    p.y <- out$p_y
  }
  return(p.y)
}

PGamma <- function(y, loc, scale, shape, c.indicator = F) {
  # Finds P(Y <= y), where Y is gamma distributed with 
  ## scale/shape parameters scale/shape.
  #
  # Args: 
  #   y: scalar.
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   P(Y <= y).
  # Error handling:
  CheckNumeric(y)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    p.y <- pgamma(y, shape = shape, scale = scale)
  } else {
    out <- .C("PGamma",
              y = as.double(y), 
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              p_y = as.double(0)
    )
    p.y <- out$p_y
  }
  return(p.y)
}

QNorm <- function(tau, loc, scale, shape = 0, c.indicator = F) {
  # Finds QY(tau), where Y is normally distributed with 
  ## location scale parameters loc/scale.
  #
  # Args: 
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value q.y, where P(Y <= q.y) = tau. 
  # Error handling:
  CheckUnitInterval(tau)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    q.tau <- qnorm(tau, loc, sd = scale)
  } else {
    out <- .C("QNorm",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              q_tau = as.double(0)
    )
    q.tau <- out$q_tau
  }
  return(q.tau)
}

QT <- function(tau, loc, scale, shape, c.indicator = F) {
  # Finds QY(tau), where Y is student's t distributed with 
  ## location/scale/df parameters loc/scale/shape.
  #
  # Args: 
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value q.y, where P(Y <= q.y) = tau. 
  # Error handling:
  CheckUnitInterval(tau)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    q.tau <- scale * qt(tau, df = shape) + loc
  } else {
    out <- .C("QT",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              q_tau = as.double(0)
    )
    q.tau <- out$q_tau
  }
  return(q.tau)
}

QLogistic <- function(tau, loc, scale, shape = 0, c.indicator = F) {
  # Finds QY(tau), where Y is logistically distributed with 
  ## location scale parameters loc/scale.
  #
  # Args: 
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: (ignored).
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value q.y, where P(Y <= q.y) = tau.
  # Error handling:
  CheckUnitInterval(tau, length.one = F)
  CheckNumeric(loc, length.one = F)
  CheckNumeric(scale, length.one = F, positive = T)
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
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    q.tau <- qlogis(tau, loc, scale)
  } else {
    if(n != 1) {
      stop(message = "C version of this function assumes scalar inputs.\n")
    }
    out <- .C("QLogistic",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(0),
              q_tau = as.double(0)
    )
    q.tau <- out$q_tau
  }
  return(q.tau)
}

QAsymmetricLaplace <- function(tau, loc, scale, shape, c.indicator = F) {
  # Finds QY(tau), where Y is asymmetric Lapace distributed with 
  ## location/scale/shape parameters loc/scale/shape.
  #
  # Args: 
  #   c.indicator: an indicator if the computation should be performed in C.
  #   tau: quantile level in (0, 1).
  #   loc: location parameter.
  #   scale: positive scale parameter.
  #   shape: shape parameter in (0, 1).
  # Returns:
  #   The value q.y, where P(Y <= q.y) = tau. 
  # Error handling:
  CheckUnitInterval(tau)
  CheckNumeric(loc)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    q.tau <- qalap(tau, loc, scale, tau = shape)
  } else {
    out <- .C("QAsymmetricLaplace",
              tau = as.double(tau),
              loc  = as.double(loc),
              scale  = as.double(scale),
              shape  = as.double(shape),
              q_tau = as.double(0)
    )
    q.tau <- out$q_tau
  }
  return(q.tau)
}

QWeibull <- function(tau, loc, scale, shape, c.indicator = F) {
  # Finds QY(tau), where Y is Weibull distributed with 
  ## scale/shape parameters scale/shape.
  #
  # Args: 
  #   tau: quantile level in (0, 1).
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value q.y, where P(Y <= q.y) = tau.
  # Error handling:
  CheckUnitInterval(tau)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    q.tau <- qweibull(tau, shape, scale)
  } else {
    out <- .C("QWeibull",
              tau = as.double(tau),
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              q_tau = as.double(0)
    )
    q.tau <- out$q_tau
  }
  return(q.tau)
}

QGamma <- function(tau, loc, scale, shape, c.indicator = F) {
  # Finds QY(tau), where Y is gamma distributed with 
  ## scale/shape parameters scale/shape. 
  #
  # Args: 
  #   tau: quantile level in (0, 1).
  #   loc: (ignored).
  #   scale: positive scale parameter.
  #   shape: positive shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   The value q.y, where P(Y <= q.y) = tau. 
  # Error handling:
  CheckUnitInterval(tau)
  CheckNumeric(scale, positive = T)
  CheckNumeric(shape, positive = T)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    q.tau <- qgamma(tau, shape = shape, scale = scale)
  } else {
    out <- .C("QGamma",
              tau = as.double(tau), 
              loc  = as.double(0),
              scale  = as.double(scale),
              shape  = as.double(shape),
              q_tau = as.double(0)
    )
    q.tau <- out$q_tau
  }
  return(q.tau)
}

TruncatedGaussianMoments <- function(mu, sigma, a, b, type, c.indicator = F) {
  # Finds mean and variance for a truncated Gaussian random variable.
  # Args: 
  #   mu: mean.
  #   sigma: standard deviation.
  #   a: left-censoring point (can be -Inf).
  #   b: right-censoring point (can be Inf).
  #   type: an indicator for the type of censoring. 
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   mn:  means of the distribution.
  #   vr:  variances of the distribution.
  # Error handling
  
  CheckNumeric(mu)
  CheckNumeric(sigma, positive = T)
  CheckNumeric(a)
  CheckNumeric(b)
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
  
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    if(type == "left") {
      a.centered <- (a - mu) / sigma
      lambda <- dnorm(a.centered) / (1 - pnorm(a.centered))
      delta <- lambda * (lambda - a.centered)
      mn <- mu + sigma * lambda 
      vr <- sigma2 * (1 - delta)      
    } else if (type == "right") {
      b.centered <- (b - mu) / sigma
      b.centered.ratio <- dnorm(b.centered) / pnorm(b.centered)
      mn <- mu - sigma * b.centered.ratio
      vr <- sigma2 * (1 - b.centered * b.centered.ratio - b.centered.ratio ^ 2)      
    } else {
      a.centered <- (a - mu) / sigma
      b.centered <- (b - mu) / sigma    
      ratio <- (dnorm(a.centered) - dnorm(b.centered)) / (pnorm(b.centered) - pnorm(a.centered))      
      mn <- mu + sigma * ratio
      vr <- sigma2 * (1 + 
                        (a.centered * dnorm(a.centered) - b.centered * dnorm(b.centered)) / 
                        (pnorm(b.centered) - pnorm(a.centered)) - 
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
    out <- .C("TruncatedGaussianMoments",
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

RandomNormal <- function(n.samples, mu, sigma, c.indicator = F) {
  # Generates a multivariate normal random variable.
  # The sample will differ depending on c.indicator, but both versions are valid.
  #
  # Args: 
  #   n.samples: number of samples.
  #   mu: the mean parameter.
  #   sigma: the covariance parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   y: multivariate normal realization.
  #
  # Error handling:
  CheckInteger(n.samples, positive = T)
  CheckNumeric(mu, length.one = F)
  CheckInvertibleSymmetricMatrix(sigma)
  dim.y <- length(mu)  
  if(nrow(sigma) != dim.y) {
    stop(message = "number of rows of sigma and length of mu must be equal.\n")
  }
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    y <- rmvnorm(n.samples, mean = mu, sigma = sigma)
  } else {
    out <- .C("RandomNormal",
              n_samples   = as.integer(n.samples),
              dim_y       = as.integer(dim.y),            
              mu          = as.double(mu),            
              sigma       = as.double(sigma),            
              y           = as.double(rep(0,  n.samples * dim.y))
    )    
    y <- array(out$y, dim = c(n.samples, dim.y))
  }
  return(y)
}

RandomTruncatedNormal <- function(n.samples, mu, sigma, a, b, c.indicator = F){
  # Generates a truncated normal random variable.
  #
  # Args: 
  #   n.samples: number of samples.
  #   mu: the mean parameter.
  #   sigma: the scale parameter.
  #   a: the left truncation point.  Can be -Inf.
  #   b: the right truncation point.  Can be Inf.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   y: truncated normal realization.
  #
  # Error handling:
  CheckInteger(n.samples, positive = T)
  CheckNumeric(mu, length.one = F)
  if(length(mu) == 1) {
    mu <- rep(mu, n.samples)
  }
  if(length(mu) != n.samples) {
    stop(message = "length of mu must be n.samples.\n")
  }
  
  CheckNumeric(sigma, positive = T, length.one = F)
  if(length(sigma) == 1) {
    sigma <- rep(sigma, n.samples)
  }
  if(length(sigma) != n.samples) {
    stop(message = "length of sigma must be n.samples.\n")
  }
  CheckNumeric(a, length.one = F)
  if(length(a) == 1) {
    a <- rep(a, n.samples)
  }
  if(length(a) != n.samples) {
    stop(message = "length of a must be n.samples.\n")
  }
  CheckNumeric(b, length.one = F)
  if(length(b) == 1) {
    b <- rep(b, n.samples)
  }
  if(length(b) != n.samples) {
    stop(message = "length of b must be n.samples.\n")
  }
  if(any(b < a)) {
    stop(message = "b must be greater than or equal to a for all elements.\n")
  }
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    u <- runif(n.samples)  
    pnorm.alpha <- ifelse(a == -Inf, 0, pnorm((a - mu) / sigma))
    pnorm.beta <- ifelse(b == Inf, 1, pnorm((b - mu) / sigma))
    y <- qnorm(pnorm.alpha + u * (pnorm.beta - pnorm.alpha)) * sigma + mu 
  } else {
    a <- ifelse(a == -Inf, -99, a)
    b <- ifelse(b ==  Inf, -99, b)    
    
    out <- .C("RandomTruncatedNormal",
              n_samples   = as.integer(n.samples),
              mu          = as.double(mu),
              sigma       = as.double(sigma),
              a           = as.double(a),          
              b           = as.double(b),
              y           = as.double(rep(0, n.samples))
    )            
    y <- matrix(out$y, ncol = 1)
  }
  return(y)
}

RandomInverseWishart <- function(df, sai, c.indicator = F) {
  # Generates an inverse Wishart random variable.
  # The sample will differ depending on c.indicator, but both versions are valid.
  #
  # Args: 
  #   df: the degrees of freedom.
  #   sai: the scale matrix.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   w: inverse Wishart realization.
  #
  # Error handling:
  
  CheckNumeric(df, positive = T)
  CheckInvertibleSymmetricMatrix(sai)
  CheckBinary(c.indicator)
  
  if(!c.indicator) {
    w <- riwish(df, sai)
  } else {
    n.row <- nrow(sai)
    out <- .C("RandomInverseWishart",
              n_row   = as.integer(n.row),
              df      = as.double(df),          
              sai     = as.double(sai),
              w       = as.double(rep(0, n.row * n.row))
    )            
    w <- array(out$w, dim = c(n.row, n.row))
  }
  return(w)
}

RandomDirichlet <- function(n.samples = 1, probs, c.indicator = F) {
  # Generates n.samples draws from a Dirichlet(probs) distribution.  
  #   
  # Args:
  #   n.samples: the number of samples.
  #   probs: vector of probabilities. 
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : n.samples x length(probs) array of Dirichlet draws. 
  # Error handling:
  CheckInteger(n.samples, positive = T)
  CheckUnitInterval(probs, length.one = F)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    x <- rdirichlet(n.samples, probs)
  } else {
    n.cells <- length(probs)    
    out <- .C("RandomDirichlet",
              n_samples       = as.integer(n.samples),            
              n_cells         = as.integer(n.cells),            
              probs           = as.double(probs),
              x               = as.double(rep(0, n.samples * n.cells))
    )
    x <- array(out$x, dim = c(n.samples, n.cells))
  } 
  return(x)
}

RandomMultinomial <- function(n.samples = 1, probs, c.indicator = F) {
  # Generates n.samples draws from a multinomial(probs) distribution.  
  #   
  # Args:
  #   n.samples: the number of samples.
  #   probs: vector of probabilities. 
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : length(probs) vector of multinomial draws. 
  # Error handling:
  CheckInteger(n.samples, positive = T)
  CheckUnitInterval(probs, length.one = F)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    x <- rmultinom(n = 1, size = n.samples, prob = probs)
  } else {
    n.cells <- length(probs)    
    out <- .C("RandomMultinomial",
              n_samples         = as.integer(n.samples),            
              n_cells           = as.integer(n.cells),
              probs             = as.double(probs),
              x                 = as.integer(rep(0, n.cells))            
    )
    x <- matrix(out$x, ncol = 1)
  }
  return(x)
}

RandomUniform <- function(a = 0, b = 1, c.indicator = F) {
  # Generates a draw from a unif(a, b) distribution.  
  #   
  # Args:
  #   a: lower limit in the support.
  #   b: upper limit in the support.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : unif(a, b) variable. 
  # Error handling:
  CheckNumeric(a)
  CheckNumeric(b)
  if(a > b) {
    stop(message = "a must be less than or equal to b.\n")
  }
  CheckBinary(c.indicator)
  if(!c.indicator) {
    x <- runif(1, a, b)
  } else {
    out <- .C("RandomUniform",
              a = as.double(a),
              b = as.double(b),
              x = as.double(0)            
    )
    x <- out$x
  }
  return(x)
}

RandomBeta <- function(shape1, shape2, c.indicator = F) {
  # Generates a draw from a beta(shape1, shape2) distribution.  
  #   
  # Args:
  #   shape1: the first shape parameter.
  #   shape2: the second shape parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : beta random variable. 
  # Error handling:
  CheckNumeric(shape1, positive = T)
  CheckNumeric(shape2, positive = T)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    x <- rbeta(n = 1, shape1, shape2)
  } else {
    out <- .C("RandomBeta",
              shape1 = as.double(shape1),
              shape2 = as.double(shape2),
              x      = as.double(0)            
    )
    x <- out$x
  }
  return(x)
}

RandomGamma <- function(shape, scale, c.indicator = F) {
  # Generates a draw from a gamma(shape, scale) distribution.  
  #   
  # Args:
  #   shape: the shape parameter.
  #   scale: the scale parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : gamma random variable. 
  # Error handling:
  CheckNumeric(shape, positive = T)
  CheckNumeric(scale, positive = T)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    x <- rgamma(n = 1, shape = shape, scale = scale)
  } else {
    n.cells <- length(probs)    
    out <- .C("RandomGamma",
              shape = as.double(shape),
              scale = as.double(scale),
              x     = as.double(0)            
    )
    x <- out$x
  }
  return(x)
}

RandomExponential <- function(scale, c.indicator = F) {
  # Generates a draw from a gamma(shape, scale) distribution.  
  #   
  # Args:
  #   scale: the scale parameter.
  #   c.indicator: an indicator if the computation should be performed in C.
  # Returns:
  #   x : exponential random variable. 
  # Error handling:
  CheckNumeric(scale, positive = T)
  CheckBinary(c.indicator)
  if(!c.indicator) {
    x <- rexp(n = 1, rate = 1 / scale)
  } else {
    out <- .C("RandomExponential",
              scale = as.double(scale),
              x     = as.double(0)            
    )
    x <- out$x
  }
  return(x)
}

UpdateKroneckerInverseWishart <- function(nu.0, sigma.0, e, omega, c.indicator = F) {
  # Performs a Gibbs update of the posterior inverse Wishart random variable sigma,
  # where e|omega, sigma ~ N(0, Kronecker(omega, sigma)) and sigma ~ IW(nu.0, sigma.0).
  #
  # Args:
  #   nrow.omega: the number of rows of omega.
  #   nrow.sigma: the number of rows of sigma.
  #   nu.0: the prior scale for sigma.
  #   sigma.0: the prior location for sigma.
  #   e: a multivariate normal random variable that has mean 0
  #     and covariance Kronecker(omega ^ -1, sigma)).
  #   omega: a valid precision matrix.
  #
  # Returns:
  #   sigma: the updated covariance.
  # Error Handling:
  CheckNumeric(nu.0, positive = T)
  CheckInvertibleSymmetricMatrix(sigma.0)
  CheckNumeric(e, length.one = F)  
  if(!isSymmetric(omega)) {
    stop(message = "omega should be a symmetric matrix.\n")
  }
  
  nrow.omega <- nrow(omega)
  nrow.sigma <- nrow(sigma.0)
  
  if(length(e) != (nrow.omega * nrow.sigma)) {
    stop(message = "e must be of length nrow(omega) * nrow(sigma.0).\n")
  }
  CheckBinary(c.indicator)
  if(!c.indicator) {
    s <- sigma.0
    for(i in 1:nrow.omega) {
      for(j in 1:nrow.omega) {
        this.e1 <- e[((i - 1) * nrow.sigma + 1) : (i * nrow.sigma)]
        this.e2 <- e[((j - 1) * nrow.sigma + 1) : (j * nrow.sigma)]
        s <- s + omega[i,j] * (this.e2 %*% t(this.e1)) 
      }
    }
    sigma <- riwish(nrow.omega + nu.0, s)
  } else {
    out <- .C("UpdateKroneckerInverseWishart",
              nrow_omega      = as.integer(nrow.omega),
              nrow_sigma      = as.integer(nrow.sigma),
              nu_0            = as.double(nu.0),
              sigma_0         = as.double(sigma.0),
              e               = as.double(e),
              omega           = as.double(omega),
              sigma           = as.double(rep(0, nrow.sigma * nrow.sigma))
              
    )
    sigma <- array(out$sigma, dim = c(nrow.sigma, nrow.sigma))
  }
  return(sigma)  
}





