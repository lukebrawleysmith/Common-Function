#setwd("C:/Users/smithinger/Desktop/bsquare_gjh/New")
#setwd("H:/My Documents/BSquareNew")
setwd("E:/New")
source("CommonFunctions.R")
source("PriorFunctions.R")
source("DataFunctions.R")
dyn.load("BSquare06042015.dll")

#check CheckCharacter, CheckBinary, CheckInteger, CheckNumeric, CheckUnitInverval,
# CheckLocs, CheckInteriorKnots     
{
  luke <- "centered"
  CheckCharacter(luke)
#  CheckCharacter(pi)
  a <- 3
#  CheckCharacter(a)

  CheckBinary(0)
  CheckBinary(1)
#  CheckBinary(c(0, 1))
  CheckBinary(c(0, 1), length.one = F)
  #CheckBinary(c(0, 1, NA), length.one = F)
  CheckBinary(c(0, 1, NA), length.one = F, no.missing = F)

  #CheckInteger(pi)
  CheckInteger(-3)
  #CheckInteger(c(-3, 0, 2), length.one = T)
  CheckInteger(c(-3, 0, 2), length.one = F)
  #CheckInteger(c(NA, 0, 2), length.one = F, no.missing = T)
  CheckInteger(c(NA, 0, 2), length.one = F, no.missing = F)
  #CheckInteger(v = c(-1, 0, 2), length.one = F, nonnegative = T)
  CheckInteger(c(-1, 0, 2), length.one = F, nonnegative = F)
  #CheckInteger(v = c(0, 2), length.one = F, positive = T)
  CheckInteger(c(0, 2), length.one = F, positive = F)
  
  CheckNumeric(pi)
  CheckNumeric(Inf)
  #CheckNumeric("merissa")
  CheckNumeric(-3)
  #CheckNumeric(c(-3, 0, 2), length.one = T)
  CheckNumeric(c(-3, 0, 2), length.one = F)
  #CheckNumeric(c(NA, 0, 2), length.one = F, no.missing = T)
  CheckNumeric(c(NA, 0, 2), length.one = F, no.missing = F)
  #CheckNumeric(v = c(-1, 0, 2), length.one = F, nonnegative = T)
  CheckNumeric(c(-1, 0, 2), length.one = F, nonnegative = F)
  #CheckNumeric(v = c(0, 2), length.one = F, positive = T)
  CheckNumeric(c(0, 2), length.one = F, positive = F)
  
  #CheckUnitInterval(runif(100))
  CheckUnitInterval(runif(100), length.one = F)
  #CheckUnitInterval(1 + runif(100), length.one = F)
  #CheckUnitInterval(c(NA, runif(100)), length.one = F, no.missing = T)
  CheckUnitInterval(c(NA, runif(100)), length.one = F, no.missing = F)
  
  locs <- matrix(runif(12), 6, 2)
  CheckLocs(locs)
  #CheckLocs(matrix(locs, 4, 3))
  locs[1,1] <- locs[2,1]
  CheckLocs(locs)
  locs[1,] <- locs[2,] 
  #CheckLocs(locs)
  
  set.seed(11212015)
  knots <- runif(5)
#  CheckInteriorKnots(knots)
  knots <- sort(knots)
  CheckInteriorKnots(knots)
  spline.knots <- MakeSplineKnots(knots, spline.df = 3)
  CheckCubicSplineKnots(spline.knots)  
}

#check SetEqual, Plot, Mean, Variance, MakeStickBreakingWeights, FindInterval
{
  
  uuu <- runif(100)
  max(abs(uuu - SetEqual(uuu)))
  max(abs(uuu - SetEqual(uuu, c.indicator = T)))
  
  u <- runif(100)
  v <- runif(100)
  
  Plot(u, v)
  
  x <- sort(runif(10))
  
  mean(x) - Mean(x, c.indicator = T)
  mean(x) - Mean(x, c.indicator = F)
  
  Variance(x[1], c.indicator = F)
  var(x) - Variance(x, c.indicator = F)
  var(x) - Variance(x, c.indicator = T)
  
  #Brian's MakeStickBreakingWeights function
  makeprobs<-function(v){ 
    N <- length(v)
    probs<-v
    probs[2:N]<-probs[2:N]*cumprod(1-v[2:N-1])
    probs}
  
  v <- x
  v[length(v)] <- 1
  max(abs(makeprobs(v) - MakeStickBreakingWeights(v, c.indicator = F)))
  max(abs(makeprobs(v) - MakeStickBreakingWeights(v, c.indicator = T)))
  
  tau <- x[1] - 2
  #(FindInterval(tau, x, c.indicator = F) - 1) - FindInterval(tau, x, c.indicator = T)
  tau <- x[1]
  (FindInterval(tau, x, c.indicator = F) - 1) - FindInterval(tau, x, c.indicator = T)
  tau <- x[1] + 10 ^ -5
  (FindInterval(tau, x, c.indicator = F) - 1) - FindInterval(tau, x, c.indicator = T)
  tau <- 0.5 * (min(x) + max(x))
  (FindInterval(tau, x, c.indicator = F) - 1) - FindInterval(tau, x, c.indicator = T)  
  tau <- x[length(x)]  
  (FindInterval(tau, x, c.indicator = F) - 1) - FindInterval(tau, x, c.indicator = T)  
  tau <- 20
  (FindInterval(tau, x, c.indicator = F) - 1) - FindInterval(tau, x, c.indicator = T)  
}

#check MatrixMultiply, Kronecker
{
  nrow.a <- 10
  ncol.a <- 5
  nrow.b <- 5
  ncol.b <- 3
  
  a <- matrix(runif(nrow.a * ncol.a), nrow.a, ncol.a)
  b <- matrix(runif(nrow.b * ncol.b), nrow.b, ncol.b)
  c <- a %*% b
  c1 <- MatrixMultiply(a, b, c.indicator = T)
  c2 <- MatrixMultiply(a, b, c.indicator = F)
  max(abs(c - c1))
  max(abs(c - c2))
  
  a <- t(a)
  c <- t(a) %*% b
  c1 <- MatrixMultiply(a, b, trans.a = T, c.indicator = T)
  c2 <- MatrixMultiply(a, b, trans.a = T, c.indicator = F)
  max(abs(c - c1))
  max(abs(c - c2))
  
  a <- t(a)
  b <- t(b)
  c <- a %*% t(b)
  c1 <- MatrixMultiply(a, b, trans.a = F, trans.b = T, c.indicator = T)
  c2 <- MatrixMultiply(a, b, trans.a = F, trans.b = T, c.indicator = F)
  max(abs(c - c1))
  max(abs(c - c2))

  a <- t(a)
  c <- t(a) %*% t(b)
  c1 <- MatrixMultiply(a, b, trans.a = T, trans.b = T, c.indicator = T)
  c2 <- MatrixMultiply(a, b, trans.a = T, trans.b = T, c.indicator = F)
  max(abs(c - c1))
  max(abs(c - c2))
  
  c <- kronecker(a, b)
  c1 <- Kronecker(a, b, c.indicator = T)
  c2 <- Kronecker(a, b, c.indicator = F)
  max(abs(c - c1))
  max(abs(c - c2))
}

#check QuadraticForm, QuadraticFormKronecker, CholMatrix, 
## InvertSymmetricMatrix, CheckInvertibleSymmetricMatrix
{
  n.responses <- 10
  v <- (n.responses) + 5
  s <- (v - n.responses - 1) * diag(n.responses)
  cov.responses <- riwish(v, s)
  
  x <- runif(n.responses)
  truth <- t(x) %*% cov.responses %*% x
  QuadraticForm(x, cov.responses, c.indicator = F) - truth
  QuadraticForm(x, cov.responses, c.indicator = T) - truth
  
  nrow.d <- 3
  d.mat <- riwish(v = nrow.d, diag(nrow.d))
  nrow.e <- 4
  e.mat <- riwish(v = nrow.e, diag(nrow.e))
  x <- runif(nrow.d * nrow.e)
  
  x <- 1:(nrow.d * nrow.e)
  
  qf <- t(x) %*% kronecker(d.mat, e.mat) %*% x
  qf2 <- QuadraticFormKronecker(d.mat, e.mat, x, c.indicator = F)
  qf3 <- QuadraticFormKronecker(d.mat, e.mat, x, c.indicator = T)
  qf - qf2
  qf - qf3
  
  truth <- t(chol(cov.responses))
  out <- CholMatrix(cov.responses, c.indicator = F)
  max(abs(out$v.chol - t(chol(cov.responses))))
  out$log.det - log(det(cov.responses))
  out <- CholMatrix(cov.responses, c.indicator = T)
  max(abs(out$v.chol - t(chol(cov.responses))))
  out$log.det - log(det(cov.responses))
  
  out <- InvertSymmetricMatrix(cov.responses, c.indicator = F)
  prec.responses <- out$v.inv
  log.det <- out$log.det
  max(abs(diag(n.responses) - prec.responses %*% cov.responses))
  log(det(solve(cov.responses))) - log.det  
  
  out <- InvertSymmetricMatrix(cov.responses, c.indicator = T)
  prec.responses.2 <- out$v.inv
  log.det.2 <- out$log.det
  max(abs(prec.responses - prec.responses.2))
  log.det - log.det.2
  
  mat.1 <- matrix(1, 2, 2)
  CheckInvertibleSymmetricMatrix(prec.responses.2)
  #CheckInvertibleSymmetricMatrix(mat.1)
}

#check MakeAbsoluteDistance, MakeAutoregressivePrecision
{
  n.timepoints <- 5
  timepoints <- 1:n.timepoints
  
  d <- matrix(NA, n.timepoints, n.timepoints)
  for(i in 1:n.timepoints) {
    for(j in 1:n.timepoints) {
      d[i, j] <- abs(timepoints[i] - timepoints[j])
    }
  }
  d2 <- MakeAbsoluteDistance(timepoints)
  d2.c <- MakeAbsoluteDistance(timepoints, c.indicator = T)
  max(abs(d - d2))  
  max(abs(d - d2.c))
  
  rho <- 0.5
  rho <- runif(1)
  
  p1 <- solve(rho ^ d)
  out2 <- MakeAutoregressivePrecision(rho, d, c.indicator = F) 
  p2 <- out2$prec.ar
  log.det2 <- out2$log.det
  
  out3 <- MakeAutoregressivePrecision(rho, d, c.indicator = T)
  p3 <- out3$prec.ar
  log.det3 <- out3$log.det
  
  max(abs(p1 - p2))
  max(abs(p1 - p3))
  
  log(det(p1)) - log.det2
  log(det(p1)) - log.det3
  
}

#check CheckLocs, MakeSpatialDistance, MakeSpatialPrecision
{
  n.locs <- 5
  locs <- matrix(runif(n.locs * 2), ncol = 2)
  
  CheckLocs(locs)
  dummy.locs <- locs
  dummy.locs[1, 1] <- dummy.locs[1 ,2] 
  CheckLocs(dummy.locs)
  dummy.locs[1, ] <- dummy.locs[2, ] 
  #CheckLocs(dummy.locs)
  
  d <- matrix(NA, n.locs, n.locs)
  for(i in 1:n.locs) {
    for(j in 1:n.locs) {
      d[i, j] <- sqrt((locs[i, 1] - locs[j, 1]) ^ 2 +  (locs[i, 2] - locs[j, 2])  ^ 2)
    }
  }
  d2 <- MakeSpatialDistance(locs)
  max(abs(d - d2))  
  
  phi <- 0.5
  phi <- runif(1)
  
  p1 <- solve(exp(-d / phi))
  out2 <- MakeSpatialPrecision(phi, d, c.indicator = F) 
  p2 <- out2$prec.sp
  log.det2 <- out2$log.det
  
  out3 <- MakeSpatialPrecision(phi, d, c.indicator = T)
  p3 <- out3$prec.sp
  log.det3 <- out3$log.det
  
  max(abs(p1 - p2))
  max(abs(p1 - p3))
  
  log(det(p1)) - log.det2
  log(det(p1)) - log.det3
}

#check MakeSplineKnots, MSpline, ISpline, CubicMSpline, CubicISpline,
# MakeSplineCoefficients, CubicMSpline2, CubicISpline2
{
  set.seed(11212015)
  spline.df <- 3
  interior.knots <- c(0.3, 0.5, 0.6)
  truth <- c(rep(0, spline.df), interior.knots, rep(1, spline.df + 2))
  spline.knots1 <- MakeSplineKnots(interior.knots, spline.df, c.indicator = F)
  spline.knots2 <- MakeSplineKnots(interior.knots, spline.df, c.indicator = T)
  max(abs(spline.knots1 - truth))
  max(abs(spline.knots2 - truth))  
  
  tau <- 0.05
  
  spline.df <- 1
  n.splines <- spline.df + length(interior.knots)
  spline.knots <- MakeSplineKnots(interior.knots, spline.df)
  findInterval(tau, MakeSplineKnots(interior.knots, spline.df))
  sapply(1:n.splines, FUN = MSpline, tau = tau, spline.df = spline.df, 
         spline.knots = spline.knots)
  
  spline.df <- 3
  n.splines <- spline.df + length(interior.knots)
  spline.knots <- MakeSplineKnots(interior.knots, spline.df)
  findInterval(tau, MakeSplineKnots(interior.knots, spline.df))
  sapply(1:n.splines, FUN = MSpline, tau = tau, spline.df = spline.df, 
         spline.knots = spline.knots)
  
  #check CubicMSpline
  #  Jingwen's mspline function
  # I modified this function so knots.inter was an argument
  D.I<-function(x, knots.inter)
  {
    boundary.knot<-c(0,1)
    origin<-floor(boundary.knot[1])      #lower level 
    end<-ceiling(boundary.knot[2])       #upper level
    
    df<-3
    n.knots<-length(knots.inter)+df
    knots.all<-c(rep(0,df),knots.inter,rep(1,df))#the knots sequents with boundary
    #all knots : n.knots+1
    nknots<-n.knots
    if( x!=boundary.knot[2])
    {
      it<-findInterval(x,knots.all)  
      I.temp<-rep(0,nknots)  #find interval for x
      #define the basis functions
      I.temp[min((it+1),nknots):nknots]<-0
      I.temp[1:max(1,(it-3))]<-0
      I.temp[it]<-3*(x-knots.all[it])^2/((knots.all[it+3]-knots.all[it])*(knots.all[it+2]-knots.all[it])*(knots.all[it+1]-knots.all[it]))
      I3<- -3*(knots.all[it+2]-x)^2/((knots.all[it+2]-knots.all[it-1])*(knots.all[it+2]-knots.all[it])*(knots.all[it+1]-knots.all[it]))-3*(x-knots.all[it-1])^2/((knots.all[it+2]-knots.all[it-1])*(knots.all[it+1]-knots.all[it-1])*(knots.all[it+1]-knots.all[it]))
      I1<-  3/(knots.all[it+1]-knots.all[it])
      I.temp[it-1]<-I1+I3
      I.temp[it-2]<-3*(knots.all[it+1]-x)^2/((knots.all[it+1]-knots.all[it-2])*(knots.all[it+1]-knots.all[it-1])*(knots.all[it+1]-knots.all[it]))
    }
    if(x==boundary.knot[2]) {I.temp<-c(rep(0,(nknots-1)),3/(x-knots.all[nknots]))}
    return(I.temp)
  }
  
  #  Jingwen's ispline function: Integrated I-spline#
  # I modified this function so knots.inter was an argument
  I.basis<-function(x, knots.inter)
  {
    
    boundary.knot<-c(0,1)
    origin<-floor(boundary.knot[1])      #lower level 
    end<-ceiling(boundary.knot[2])       #upper level
    
    df<-3
    n.knots<-length(knots.inter)+df
    knots.all<-c(rep(0,df),knots.inter,rep(1,df))#the knots sequents with boundary
    #all knots : n.knots+1
    nknots<-n.knots
    if( x!=boundary.knot[2])
    {
      it<-findInterval(x,knots.all)  
      I.temp<-rep(0,nknots)  #find interval for x
      #define the basis functions
      I.temp[min((it+1),nknots):nknots]<-0
      I.temp[1:max(1,(it-3))]<-1
      I.temp[it]<-(x-knots.all[it])^3/((knots.all[it+3]-knots.all[it])*(knots.all[it+2]-knots.all[it])*(knots.all[it+1]-knots.all[it]))
      I3<- (knots.all[it+2]-x)^3/((knots.all[it+2]-knots.all[it-1])*(knots.all[it+2]-knots.all[it])*(knots.all[it+1]-knots.all[it]))-(x-knots.all[it-1])^3/((knots.all[it+2]-knots.all[it-1])*(knots.all[it+1]-knots.all[it-1])*(knots.all[it+1]-knots.all[it]))
      I2<-  (x-knots.all[it])^2/((knots.all[it+2]-knots.all[it])*(knots.all[it+1]-knots.all[it]))-(knots.all[it+2]-x)^2/((knots.all[it+2]-knots.all[it])*(knots.all[it+1]-knots.all[it]))
      I1<-  (x-knots.all[it-1])/(knots.all[it+1]-knots.all[it])
      I.temp[it-1]<-I1+I2+I3
      I.temp[it-2]<-1-(knots.all[it+1]-x)^3/((knots.all[it+1]-knots.all[it-2])*(knots.all[it+1]-knots.all[it-1])*(knots.all[it+1]-knots.all[it]))
    }
    if(x==boundary.knot[2]) {I.temp<-rep(1,nknots)}
    return(I.temp)
  }
  
  c1 <- MakeCubicSplineCoefficients(interior.knots, c.indicator = F)
  c2 <- MakeCubicSplineCoefficients(interior.knots, c.indicator = T)
  max(abs(c1 - c2))
  
  
  tau <- runif(100)
  tau <- c(tau, interior.knots)
  tau <- sort(tau)
  
  v.jw <- v <- vc <- v2 <- v2c <- v3 <- v3c <- matrix(NA, length(tau), 3 + length(interior.knots))
  w.jw <- w <- wc <- w2 <- w2c <- w3 <- w3c <- matrix(NA, length(tau), 3 + length(interior.knots))
  
  for(i in 1:length(tau)) {
    v.jw[i, ] <- D.I(tau[i], interior.knots)
    w.jw[i, ] <- I.basis(tau[i], interior.knots)
    for(m in 1:ncol(v)) {
      v[i, m] <- CubicMSpline(tau[i], interior.knots, m, c.indicator = F)
      vc[i, m] <- CubicMSpline(tau[i], interior.knots, m, c.indicator = T)
      v2[i, m] <- MSpline(tau[i], spline.df = 3, spline.knots, m, c.indicator = F)
      v2c[i, m] <- MSpline(tau[i], spline.df = 3, spline.knots, m, c.indicator = T)      
      w[i, m] <- CubicISpline(tau[i], interior.knots, m, c.indicator = F)
      wc[i, m] <- CubicISpline(tau[i], interior.knots, m, c.indicator = T)
      w2[i, m] <- ISpline(tau = tau[i], spline.df = 3, spline.knots, m, c.indicator = F)
      w2c[i, m] <- ISpline(tau[i], spline.df = 3, spline.knots, m, c.indicator = T)
      v3[i, m] <- CubicMSpline2(tau = tau[i], spline.knots, spline.coefs = c1, m, c.indicator = F)
      v3c[i, m] <- CubicMSpline2(tau = tau[i], spline.knots, spline.coefs = c1, m, c.indicator = T)
      w3[i, m] <- CubicISpline2(tau = tau[i], spline.knots, spline.coefs = c1, m, c.indicator = F)
      w3c[i, m] <- CubicISpline2(tau = tau[i], spline.knots, spline.coefs = c1, m, c.indicator = T)
    }
  }

  max(abs(v.jw - v))
  max(abs(v.jw - vc))
  max(abs(v.jw - v2))
  max(abs(v.jw - v2c))
  max(abs(v.jw - v3))
  max(abs(v.jw - v3c))
  max(abs(w.jw - w))
  max(abs(w.jw - wc))
  max(abs(w.jw - w2))
  max(abs(w.jw - w2c))
  max(abs(w.jw - w3))
  max(abs(w.jw - w3c))
}

#check LogDNorm, LogDT, LogDLogistic, LogDAsymmetricLapace, LogDWeibull, LogDGamma,
#      PNorm, PT, PLogistic, PAsymmetricLapace, PWeibull, PGamma,
#      QNorm, QT, QLogistic, QAsymmetricLapace, QWeibull, QGamma
{
  library(VGAM)
  set.seed(11092015)
  y <- rnorm(1)
  loc <- rnorm(1)
  scale <- rgamma(1, shape = 2, scale = 2)
  shape <- rgamma(1, shape = 2, scale = 2)
  shape2 <- runif(1)
  z <- (y - loc) / scale
  tau <- runif(1)
  
  LogDNorm(y, loc, scale, shape, c.indicator = F)  - dnorm(y, loc, sd = scale, log = T)
  LogDT(y, loc, scale, shape, c.indicator = F) - (dt(z, df = shape, log = T) - log(scale))
  LogDLogistic(y, loc, scale, shape, c.indicator = F)  - dlogis(y, loc, scale, log = T)
  LogDAsymmetricLaplace(y, loc, scale, shape2, c.indicator = F)  - log(dalap(y, location = loc, scale = scale, tau = shape2))
  LogDAsymmetricLaplace(loc, y, scale, shape2, c.indicator = F)  - log(dalap(loc, location = y, scale = scale, tau = shape2))
  LogDWeibull(y, loc, scale, shape, c.indicator = F)  - dweibull(y, shape, scale, log = T)
  LogDGamma(y + loc, loc, scale, shape, c.indicator = F)  - dgamma(y + loc, shape = shape, scale = scale, log = T)
  LogDBeta(tau, shape1 = scale, shape2 = shape, c.indicator = F) - dbeta(tau, shape1 = scale, shape2 = shape, log = T)
  
  LogDNorm(y, loc, scale, shape, c.indicator = T)  - dnorm(y, loc, sd = scale, log = T)
  LogDT(y, loc, scale, shape, c.indicator = T) - (dt(z, df = shape, log = T) - log(scale))
  LogDLogistic(y, loc, scale, shape, c.indicator = T)  - dlogis(y, loc, scale, log = T)
  LogDAsymmetricLaplace(y, loc, scale, shape2, c.indicator = T)  - log(dalap(y, location = loc, scale = scale, tau = shape2))
  LogDAsymmetricLaplace(loc, y, scale, shape2, c.indicator = T)  - log(dalap(loc, location = y, scale = scale, tau = shape2))
  LogDWeibull(y, loc, scale, shape, c.indicator = T)  - dweibull(y, shape, scale, log = T)
  LogDGamma(y + loc, loc, scale, shape, c.indicator = T)  - dgamma(y + loc, shape = shape, scale = scale, log = T)
  LogDBeta(tau, shape1 = scale, shape2 = shape, c.indicator = T) - dbeta(tau, shape1 = scale, shape2 = shape, log = T)

  PNorm(y, loc, scale, shape, c.indicator = F)  - pnorm(y, loc, sd = scale)
  PT(y, loc, scale, shape, c.indicator = F)  - pt((y - loc) / scale, df = shape)
  PLogistic(y, loc, scale, shape, c.indicator = F)  - plogis(y, loc, scale)
  PAsymmetricLaplace(y, loc, scale, shape2, c.indicator = F)  - palap(y, location = loc, scale = scale, tau = shape2)
  PAsymmetricLaplace(loc, y, scale, shape2, c.indicator = F)  - palap(loc, location = y, scale = scale, tau = shape2)
  PWeibull(loc - y, loc, scale, shape, c.indicator = F)  - pweibull(loc - y, shape, scale)
  PGamma(y + loc, loc, scale, shape, c.indicator = F)  - pgamma(y + loc, shape = shape, scale = scale)
  
  PNorm(y, loc, scale, shape, c.indicator = T)  - pnorm(y, loc, sd = scale)
  PT(y, loc, scale, shape, c.indicator = T)  - pt((y - loc) / scale, df = shape)
  PLogistic(y, loc, scale, shape, c.indicator = T)  - plogis(y, loc, scale)
  PAsymmetricLaplace(y, loc, scale, shape2, c.indicator = T)  - palap(y, location = loc, scale = scale, tau = shape2)
  PAsymmetricLaplace(loc, y, scale, shape2, c.indicator = T)  - palap(loc, location = y, scale = scale, tau = shape2)
  PWeibull(loc - y, loc, scale, shape, c.indicator = T)  - pweibull(loc - y, shape, scale)
  PGamma(y + loc, loc, scale, shape, c.indicator = T)  - pgamma(y + loc, shape = shape, scale = scale)
  
  QNorm(tau, loc, scale, shape, c.indicator = F)  - qnorm(tau, loc, sd = scale)
  QT(tau, loc, scale, shape, c.indicator = F) - (scale * qt(tau, df = shape) + loc)
  QLogistic(tau, loc, scale, shape, c.indicator = F)  - qlogis(tau, loc, scale)
  QAsymmetricLaplace(tau, loc, scale, shape2, c.indicator = F)  - qalap(tau, location = loc, scale = scale, tau = shape2)
  QAsymmetricLaplace(tau, y, scale, shape2, c.indicator = F)  - qalap(tau, location = y, scale = scale, tau = shape2)
  QWeibull(tau, loc, scale, shape, c.indicator = F)  - qweibull(tau, shape, scale)
  QGamma(tau, loc, scale, shape, c.indicator = F)  - qgamma(tau, shape = shape, scale = scale)
  
  QNorm(tau, loc, scale, shape, c.indicator = T)  - qnorm(tau, loc, sd = scale)
  QT(tau, loc, scale, shape, c.indicator = T) - (scale * qt(tau, df = shape) + loc)
  QLogistic(tau, loc, scale, shape, c.indicator = T)  - qlogis(tau, loc, scale)
  QAsymmetricLaplace(tau, loc, scale, shape2, c.indicator = T)  - qalap(tau, location = loc, scale = scale, tau = shape2)
  QAsymmetricLaplace(tau, y, scale, shape2, c.indicator = T)  - qalap(tau, location = y, scale = scale, tau = shape2)
  QWeibull(tau, loc, scale, shape, c.indicator = T)  - qweibull(tau, shape, scale)
  QGamma(tau, loc, scale, shape, c.indicator = T)  - qgamma(tau, shape = shape, scale = scale)
}

#check LogDInverseWishart
{
  p <- 3
  df <- 7
  sai <- riwish(df, diag(p))
  x <- riwish(df, sai)
  
  ll <- LogDInverseWishart(df, sai, x)
  ll.c <- LogDInverseWishart(df, sai, x, c.indicator = T)
  ll - log(diwish(x, df, sai))
  ll.c - log(diwish(x, df, sai))
  
}

#check RandomNormal
{
  dim.y <- 15
  mu <- matrix(1:dim.y, ncol = 1)
  sigma <- riwish(v = 2 * dim.y , S = diag(dim.y))
  sigma.chol <- CholMatrix(sigma)$v.chol
  
  set.seed(02051980)
  y1 <- sigma.chol %*% rnorm(dim.y) + mu
  y1 <- t(y1)
  set.seed(02051980)
  y2 <- RandomNormal(n.samples = 1, mu, sigma, c.indicator = F)
  set.seed(02051980)
  y3 <- RandomNormal(n.samples = 1, mu, sigma, c.indicator = T)
  
  max(abs(y1 - y2))
  max(abs(y1 - y3))
  max(abs(y2 - y3))
  
  n.samples <- 100000

  set.seed(10062015)
  y <- matrix(NA, n.samples, dim.y)
  for(i in 1:n.samples) {
    y[i,] <- sigma.chol %*% rnorm(dim.y) + mu
  }  
  set.seed(10062015)
  y2 <- RandomNormal(n.samples, mu, sigma, c.indicator = F)
  set.seed(10062015)
  y3 <- RandomNormal(n.samples, mu, sigma, c.indicator = T)

  max(abs(y - y2)) 
  max(abs(y - y3))
  max(abs(y2 - y3))

  max(abs(apply(y, 2, mean) - mu))
  max(abs(apply(y2, 2, mean) - mu))
  max(abs(apply(y3, 2, mean) - mu))
  
  max(abs(cov(y) - sigma))
  max(abs(cov(y2) - sigma))
  max(abs(cov(y3) - sigma))
}

#check RandomTruncatedNormal, TruncatedGaussianMoments
{
  n.samples <- 10 ^ 5

  #check that standard case results in standard normal distribution
  a <- -Inf
  b <- Inf
  set.seed(11142015)
  z1 <- RandomTruncatedNormal(n.samples, 0, 1, 
                              a, b, c.indicator = F)  
  set.seed(11142015)
  z2 <- RandomTruncatedNormal(n.samples, 0, 1, 
                              a, b, c.indicator = T) 
  
  max(abs(z1 - z2))
  qqnorm(z2)
  
  #TruncatedGaussianMoments(mu = 0, sigma = 1, a, b, type = "left")
  #TruncatedGaussianMoments(mu = 0, sigma = 1, a, b, type = "right")
  #TruncatedGaussianMoments(mu = 0, sigma = 1, a, b, type = "interval")
  TruncatedGaussianMoments(mu = 0, sigma = 1, a = -2, b = 2, type = "left", c.indicator = F)
  TruncatedGaussianMoments(mu = 0, sigma = 1, a = -2, b = 2, type = "left", c.indicator = T)
  TruncatedGaussianMoments(mu = 0, sigma = 1, a = -2, b = 2, type = "right", c.indicator = F)
  TruncatedGaussianMoments(mu = 0, sigma = 1, a = -2, b = 2, type = "right", c.indicator = T)
  TruncatedGaussianMoments(mu = 0, sigma = 1, a = -2, b = 2, type = "interval", c.indicator = F)
  TruncatedGaussianMoments(mu = 0, sigma = 1, a = -2, b = 2, type = "interval", c.indicator = T)
  #add non-standard parameters
  mu <- -2
  sigma <- sqrt(2)
  #generate normals for testing against
  set.seed(11142015)
  z <- rnorm(10 ^ 5, mean = mu, sd = sigma)
  
  set.seed(11142015)
  z1 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = F)  
  set.seed(11142015)
  z2 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = T) 
  
  max(abs(z1 - z2))
  qqnorm((z2 - mu) / sigma)
  
  #left-censored data
  a <- 0
  b <- Inf
  set.seed(11142015)
  z1 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = F)
  set.seed(11142015)
  z2 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = T)
  
  max(abs(z1 - z2))
  
  these.guys <- which(z > a)
  these.guys <- intersect(these.guys, which(z < b))
  
  plot(density(z1))
  points(density(z[these.guys]), col= "green", type = "l")
  
  out <- TruncatedGaussianMoments(mu, sigma, a, b, type = "left", c.indicator = F)
  out2 <- TruncatedGaussianMoments(mu, sigma, a, b, type = "left", c.indicator = T)
  out$mn - mean(z1)
  out$v - var(z1)
  out$mn - out2$mn
  out$v - out2$v
  #interval-censored data
  a <- 0
  b <- 2
  set.seed(11142015)
  z1 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = F)
  set.seed(11142015)
  z2 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = T)
  
  max(abs(z1 - z2))
  
  these.guys <- which(z > a)
  these.guys <- intersect(these.guys, which(z < b))
  
  plot(density(z1))
  points(density(z[these.guys]), col= "green", type = "l")
  
  out <- TruncatedGaussianMoments(mu, sigma, a, b, type = "interval", c.indicator = F)
  out2 <- TruncatedGaussianMoments(mu, sigma, a, b, type = "interval", c.indicator = T)
  out$mn - mean(z1)
  out$v - var(z1)
  out$mn - out2$mn
  out$v - out2$v
  #right-censored data
  a <- -Inf
  b <- 2

  set.seed(11142015)
  z1 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = F)
  set.seed(11142015)
  z2 <- RandomTruncatedNormal(n.samples, mu, sigma, 
                              a, b, c.indicator = T)
  
  max(abs(z1 - z2))
  
  these.guys <- which(z > a)
  these.guys <- intersect(these.guys, which(z < b))
  
  plot(density(z1))
  points(density(z[these.guys]), col= "green", type = "l")
  
  out <- TruncatedGaussianMoments(mu, sigma, a, b, type = "right", c.indicator = F)
  out2 <- TruncatedGaussianMoments(mu, sigma, a, b, type = "right", c.indicator = T)
  out$mn - mean(z1)
  out$v - var(z1)
  out$mn - out2$mn
  out$v - out2$v
}

#check RandomInverseWishart, UpdateKroneckerInverseWishart
{
  nrow.sai <- 3
#  df <- 2 * nrow.sai + 1
  df <- nrow.sai + 5
  sai <- (df - nrow.sai - 1) * diag(nrow.sai)

  mean.w <- sai / (df - nrow.sai - 1)
  var.w <- matrix(NA, nrow(sai), ncol(sai))
  
  for(i in 1:nrow.sai) {
    for(j in 1:nrow.sai) {
      var.w[i, j] <- 
        ((df - nrow.sai + 1) * sai[i, j] ^ 2 + 
           (df - nrow.sai - 1) * sai[i, i] * sai[j, j]) / 
        ((df - nrow.sai) * (df - nrow.sai - 1) ^ 2 * (df - nrow.sai - 3))
    }
  }
  
  set.seed(11142015)
  v1 <- riwish(df, sai)
  set.seed(11142015)
  v2 <- RandomInverseWishart(df, sai, c.indicator = F)
  set.seed(11142015)
  v3 <- RandomInverseWishart(df, sai, c.indicator = T)

  max(abs(v1 - v2))
  max(abs(v1 - v3))

  n.iters <- 100000
  w.array <- array(NA, dim = c(n.iters, nrow.sai, nrow.sai))

  for(i in 1:n.iters) {
    w.array[i, , ] <- RandomInverseWishart(df, sai, c.indicator = T)
  }
  
  max(abs(apply(w.array, c(2, 3), mean) - mean.w))
  max(abs(apply(w.array, c(2, 3), var) - var.w))
  
  nrow.omega <- 10
  omega.df <- 20
  omega.sai <- (omega.df - nrow.omega - 1) * diag(nrow.omega)
  omega <- riwish(omega.df, omega.sai)
  omega.inv <- InvertSymmetricMatrix(omega)$v.inv
  
  nrow.sai <- 2
  df <- 2 * nrow.sai + 1
  df <- 100
  sai <- (df - nrow.sai - 1) * diag(nrow.sai)
  
  w <- RandomInverseWishart(df, sai, c.indicator = T)
  v <- kronecker(omega.inv, w)
  chol.v <- CholMatrix(v)$v.chol
  e <- chol.v %*% rnorm(nrow.omega * nrow(w))
  e <- e - mean(e)
  e.mat <- matrix(e, ncol = 2, byrow = T)

  set.seed(11142015)
  sigma1 <- UpdateKroneckerInverseWishart(df, sai, e, omega, c.indicator = F)
  set.seed(11142015)
  sigma2 <- UpdateKroneckerInverseWishart(df, sai, e, omega, c.indicator = T)
  
  max(abs(sigma1 - sigma2))

  sigma2 - w
  
  
  n.iters <- 1000
  w.array <- array(NA, dim = c(n.iters, nrow.sai, nrow.sai))
  for(i in 1:n.iters) {
    w.array[i, , ] <- UpdateKroneckerInverseWishart(df, sai, e, omega, c.indicator = T)
  }
  
  max(abs(apply(w.array, c(2, 3), mean) - w))

  #increase the dimension of omega to check for consistency
  nrow.omega <- 2000
  omega.df <- 10000
  omega.sai <- (omega.df - nrow.omega - 1) * diag(nrow.omega)
  omega <- riwish(omega.df, omega.sai)
  omega.inv <- InvertSymmetricMatrix(omega)$v.inv
  
  nrow.sai <- 2
  df <- 10
  sai <- (df - nrow.sai - 1) * diag(nrow.sai)
  
  w <- RandomInverseWishart(df, sai, c.indicator = T)
  v <- kronecker(omega.inv, w)
  chol.v <- CholMatrix(v)$v.chol
  e <- chol.v %*% rnorm(nrow.omega * nrow(w))
  e <- e - mean(e)

  set.seed(11142015)
  sigma1 <- UpdateKroneckerInverseWishart(df, sai, e, omega, c.indicator = T)
  set.seed(11142015)  
  sigma2 <- UpdateKroneckerInverseWishart(df, sai, e, omega, c.indicator = T)

  max(abs(w - sigma2))
  max(abs(sigma1 - sigma2))
}

#check RandomDirichlet, RandomMultinomial, RandomBeta, RandomGamma, RandomExponential 
{
  probs <- rdirichlet(1, runif(5))
  n.samples <- 10
  set.seed(12012015)
  x <- rdirichlet(n.samples, probs)
  set.seed(12012015)
  x2 <- RandomDirichlet(n.samples, probs, c.indicator = F)  
  set.seed(12012015)
  x3 <- RandomDirichlet(n.samples, probs, c.indicator = T)  

  max(abs(x - x2))
  max(abs(x - x3))
  
  n.samples <- 1000
  set.seed(12012015)
  x <- rmultinom(n = 1, size = n.samples, prob = probs)
  set.seed(12012015)
  x2 <- RandomMultinomial(n.samples, probs, c.indicator = F)  
  set.seed(12012015)
  x3 <- RandomMultinomial(n.samples, probs, c.indicator = T)  
  
  max(abs(x - x2))
  max(abs(x - x3))

  shape1 <- rgamma(1, shape = 1, scale = 1)
  shape2 <- rgamma(1, shape = 1, scale = 1)
  
  RandomUniform(-2, -3)
  set.seed(3)
  x1 <- RandomUniform(3, 3)
  set.seed(3)
  x2 <- RandomUniform(3, 3, c.indicator = T)
  
  x1 - x2
  
  set.seed(3)
  x1 <- RandomUniform(3, pi)
  set.seed(3)
  x2 <- RandomUniform(3, pi, c.indicator = T)
  x1 - x2
  
  
  set.seed(25)
  x1 <- rbeta(1, shape1, shape2)
  set.seed(25)
  x2 <- RandomBeta(shape1, shape2, c.indicator = F)
  set.seed(25)
  x3 <- RandomBeta(shape1, shape2, c.indicator = T)
  
  x1 - x2
  x1 - x3
  
  set.seed(10)
  x1 <- rgamma(n = 1, shape = shape1, scale = shape2)
  set.seed(10)
  x2 <- RandomGamma(shape = shape1, scale = shape2, c.indicator = F)
  set.seed(10)
  x3 <- RandomGamma(shape = shape1, scale = shape2, c.indicator = T)

  x1 - x2
  x1 - x3
  
  this.scale <- pi
  set.seed(100)
  x1 <- rexp(1, rate = 1 / this.scale)
  set.seed(100)
  x2 <- RandomExponential(scale = this.scale, c.indicator = F)
  set.seed(100)
  x3 <- RandomExponential(scale = this.scale, c.indicator = T)
  
  x1 - x2
  x1 - x3
}







