# library(sandwich)
# library(MASS)

########################### Inverting Matrices #################################

qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  mat <- crossprod(x)

  invMatrix <- tryCatch({
    chol2inv(chol(mat))
  },
  error = function(e) {
    if (grepl("leading minor of order", e$message)) {
      # If error is due to non-invertibility, issue warning and use generalized inverse
      warning("Calucating inverse of (t(X)%*%X), matrix is not positive-definite. Using generalized inverse.")
      return(ginv(mat))
    } else {
      # If it's another error, just stop and propagate the error
      stop(e)
    }
  })

  return(invMatrix)
}

########################### Weight Calculation #################################

W.fun = function(u,kernel){
  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  if (kernel.type=="Epanechnikov") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel.type=="Uniform")      w =          0.5*(abs(u)<=1)
  if (kernel.type=="Triangular")   w =   (1-abs(u))*(abs(u)<=1)
  if (kernel.type=="Gaussian")     w =   dnorm(u)
  return(w)
}

##################### Rule of Thumb Bandwidth Selection ########################

rdbw2d_rot <- function(x,kernel.type, M){

  mu2K.squared <- NA
  l2K.squared <- NA

  if (kernel.type == "Epanechnikov"){mu2K.squared <- 1/6; l2K.squared <- 4/(3 * pi)}
  if (kernel.type == "Triangular"){mu2K.squared <- 3/20; l2K.squared <- 3/(2 * pi)}
  if (kernel.type == "Uniform"){mu2K.squared <- 1/4; l2K.squared <- 1/(pi)}
  if (kernel.type == "Gaussian"){mu2K.squared <- 1; l2K.squared <- 1/(4 * pi)}

  # Data cleaning

  x <- x[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(x$x.1) & complete.cases(x$x.2)
  x <- x[na.ok,]
  N <- dim(x)[1]

  # Estimate sample variance.

  cov.matrix <- cov(x[,c("x.1", "x.2")])

  D <- 2

  trace.const <- 1/( 2^(D+2) * pi^(D/2) * det(sqrtm(cov.matrix)) ) * ( 2*sum(diag(ginv(cov.matrix, 1e-20) %*% ginv(cov.matrix, 1e-20)))
                                                                       + (sum(diag(ginv(cov.matrix, 1e-20))))^2 )

  if (is.null(M)){
    hROT <-( (D * l2K.squared) / (N * mu2K.squared * trace.const) )^(1/(4+D))
  } else{
    hROT <-( (D * l2K.squared) / (M * mu2K.squared * trace.const) )^(1/(4+D))   # Adjust for mass points.
  }

  return(hROT)
}

###### Bandwidth Selection using Global Fit for Higher Order Derivatives #######

rdbw2d_bw_v2 <- function(dat.centered, p, vec, dn, bn.1, bn.2 = NULL, vce, kernel, kernel_type, C){

  dat.centered <- dat.centered[,c("x.1", "x.2", "y", "d", "dist")]

  # Variance and coefficients for a linear combination of (p+1)-th derivatives.
  if (kernel_type == "prod"){
    w.v <- W.fun(dat.centered$x.1/c(dn), kernel) * W.fun(dat.centered$x.2/c(dn), kernel) / c(dn^2)
  }
  else{
    w.v <- W.fun(dat.centered$dist/c(dn), kernel)/c(dn^2)
  }

  ind.v <- as.logical(w.v > 0)
  eN.v <- sum(ind.v)

  ew.v <- w.v[ind.v]
  eY.v <- dat.centered$y[ind.v]

  eC.v <- C[ind.v]

  eu.v <- dat.centered[ind.v, c("x.1", "x.2")]
  eu.v$x.1 <- eu.v$x.1/dn
  eu.v$x.2 <- eu.v$x.2/dn

  if (is.null(bn.2)){
    eR.v.aug <- as.matrix(get_basis(eu.v,p+1))
    eR.v <- eR.v.aug[,1: (factorial(p+2)/(factorial(p)*2))]
    eS.v <- eR.v.aug[, (factorial(p+2)/(factorial(p)*2)+1) : (factorial(p+1+2)/(factorial(p+1)*2))]
  } else {
    eR.v.aug <- as.matrix(get_basis(eu.v,p+2))
    eR.v <- eR.v.aug[,1: (factorial(p+2)/(factorial(p)*2))]
    eS.v <- eR.v.aug[, (factorial(p+2)/(factorial(p)*2)+1) : (factorial(p+1+2)/(factorial(p+1)*2))]
    eT.v <- eR.v.aug[, (factorial(p+1+2)/(factorial(p+1)*2) + 1) : (factorial(p+2+2)/(factorial(p+2)*2))]
  }

  sqrtw_R.v <- sqrt(ew.v) * eR.v
  sqrtw_eS.v <- sqrt(ew.v) * eS.v
  sqrtw_Y.v <- sqrt(ew.v) * eY.v

  w_R.v <- ew.v * eR.v

  invG.v <- qrXXinv(sqrtw_R.v)

  vec.q <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eS.v
  vec.q <- vec.q[1,]
  vec.q <- c(rep(0, factorial(p + 2)/(factorial(p)*2)), vec.q)

  if (!is.null(bn.2)){
    sqrtw_eT.v <- sqrt(ew.v) * eT.v
    vec.t <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eT.v
    vec.t <- vec.t[1,]
    vec.t <- c(rep(0, factorial(p + 1 + 2)/(factorial(p + 1)*2)), vec.t)
  }

  invH.p <- get_invH(dn,p)

  H.p <- get_H(dn,p)

  beta.v <- invH.p %*% invG.v %*% t(sqrtw_R.v) %*% matrix(sqrtw_Y.v, ncol = 1)

  resd.v <- (eY.v - eR.v %*% (H.p %*% beta.v))[,1]

  lambda <- function(x){infl(x, invG.v)}

  if (vce=="hc0") {
    w.vce = 1
  } else if (vce=="hc1") {
    w.vce = sqrt(eN.v/(eN.v-factorial(p+2)/(factorial(p) * 2)))
  } else if (vce=="hc2") {
    hii <- apply(sqrtw_R.v, 1, lambda)
    w.vce = sqrt(1/(1-hii))
  } else if (vce=="hc3"){
    hii <- apply(sqrtw_R.v, 1, lambda)
    w.vce = 1/(1-hii)
  }

  resd.v <- resd.v * w.vce

  # sigma.v <-  t(resd.v * sqrtw_R.v) %*% (ew.v * resd.v * sqrtw_R.v) * dn^2

  sigma.v <- rd2d_vce(w_R.v, resd.v, eC.v, dn)

  V.V <- t(as.matrix(vec)) %*% t(invG.v) %*% sigma.v %*% invG.v %*% as.matrix(vec)
  V.V <- V.V[1,1]

  # Bias

  fit.ppls1 <- rd2d_lm(dat.centered, bn.1, p + 1, vce, kernel = kernel,
                       kernel_type = kernel_type, C = C, varr = TRUE)

  deriv.ppls1 <- fit.ppls1$beta
  B.B <- vec.q %*% deriv.ppls1
  B.B <- B.B[1,1]
  V.B <- matrix(vec.q, nrow = 1) %*% fit.ppls1$cov.const %*% matrix(vec.q, ncol = 1) / (bn.1^(2 + 2 * (p+1)))
  V.B <- V.B[1,1]

  Reg.v <- V.B

  Reg.b <- NA

  if (!is.null(bn.2)){
    fit.ppls2 <- rd2d_lm(dat.centered, bn.2, p + 2, vce, kernel = kernel,
                         kernel_type = kernel_type, C = C, varr = FALSE)
    deriv.ppls2 <- fit.ppls2$beta
    Reg.b <- dn * vec.t %*% deriv.ppls2
  }

  return(list("B" = B.B, "V" = V.V, "Reg.2" = Reg.b, "Reg.1" = Reg.v))
}


##################### Fitting Local Polynomial Estimators ######################

# kernel_type = "prod" or "rad"

# o = 2 means taking out all second order derivatives.

# hgrid.0 either neval or neval by 2.

# If user does not provide hgrid.1, use a single bandwidth for both treatment and control.

rd2d_fit_v2 <- function(dat, eval, deriv = NULL,o = 0, p = 1, hgrid.0, hgrid.1 = NULL,
                        kernel = "epa", kernel_type = "prod", vce = "hc1",
                        masspoints = "adjust", C = NULL, bwcheck = 50 + p + 1, unique = NULL){

  dat <- dat[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]

  sd.x1 <- sd(dat$x.1)
  sd.x2 <- sd(dat$x.2)

  neval <- dim(eval)[1]

  # Standardize: if hgrid.0 is a vector, convert it to a m by 1 data frame
  hgrid.0 <- as.data.frame(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.data.frame(hgrid.1)
  }

  if (ncol(hgrid.0) == 1){
    results <- data.frame(matrix(NA, ncol = 10, nrow = neval))
    colnames(results) <- c('ev.x.1', 'ev.x.2', 'h.0', 'h.1', 'mu.0', 'mu.1', 'se.0',
                           'se.1', 'eN.0', 'eN.1')
  } else {
    results <- data.frame(matrix(NA, ncol = 12, nrow = neval))
    colnames(results) <- c('ev.x.1', 'ev.x.2', 'h.0.x', 'h.0.y', 'h.1.x', 'h.1.y',
                           'mu.0', 'mu.1', 'se.0', 'se.1', 'eN.0', 'eN.1')
  }

  # Check for compatibility

  for (i in 1:neval){

    ev <- eval[i,]
    vec <- deriv[i,]
    h.0 <- hgrid.0[i,]
    h.1 <- h.0
    if (!is.null(hgrid.1)) h.1 <- hgrid.1[i,]

    # Center data

    dat.centered <- dat[,c("x.1", "x.2", "y", "d")]
    dat.centered$x.1 <- dat.centered$x.1 - ev$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev$x.2
    # for product kernel, standardize the covariates so that they have the same sd
    dat.centered$dist <- pmax(abs(dat.centered$x.1/sd.x1), abs(dat.centered$x.2/sd.x2)) # infinity norm
    if (kernel_type == "rad") dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2) # Euclidean norm

    # Bandwidth restriction

    if (masspoints == "adjust"){
      unique.centered <- unique
      unique.centered$x.1 <- unique.centered$x.1 - ev$x.1
      unique.centered$x.2 <- unique.centered$x.2 - ev$x.2
      unique.centered$dist <- pmax(abs(unique.centered$x.1/sd.x1), abs(unique.centered$x.2/sd.x2)) # infinity norm
      if (kernel_type == "rad") unique.centered$dist <- sqrt(unique.centered$x.1^2 + unique.centered$x.2^2) # Euclidean norm
    }

    if (!is.null(bwcheck)){

      if (masspoints == "adjust"){
        sorted.0 <- sort(unique.centered[unique.centered$d == FALSE,]$dist)
        sorted.1 <- sort(unique.centered[unique.centered$d == TRUE,]$dist)
      } else{
        sorted.0   <- sort(dat.centered[dat.centered$d == FALSE,]$dist)
        sorted.1   <- sort(dat.centered[dat.centered$d == TRUE,]$dist)
      }

      bw.min.0   <- sorted.0[bwcheck]
      bw.min.1   <- sorted.1[bwcheck]
      bw.max.0   <- sorted.0[length(sorted.0)]
      bw.max.1   <- sorted.1[length(sorted.1)]

      # convert to the original if using product kernel
      if (kernel_type == "prod"){
        multiplier <- c(sd.x1, sd.x2)
      } else {
        multiplier <- c(1,1)
      }

      bw.min.0 <- bw.min.0 * multiplier
      bw.min.1 <- bw.min.1 * multiplier
      bw.max.0 <- bw.max.0 * multiplier
      bw.max.1 <- bw.max.1 * multiplier

      if (!is.null(hgrid.1)){
        h.0     <- pmax(h.0, bw.min.0)
        h.1     <- pmax(h.1, bw.min.1)
        h.0     <- pmin(h.0, bw.max.0)
        h.1     <- pmin(h.1, bw.max.1)
      } else{
        h.0 <- pmax(h.0, bw.min.0,bw.min.1)
        h.0 <- pmin(h.0, pmax(bw.max.0, bw.max.1))
        h.1 <- h.0
      }
    }

    fit.0.p <- rd2d_lm(dat.centered[dat.centered$d == 0,], h.0, p, vce, kernel = kernel, C = C[dat.centered$d == 0],
                       varr = TRUE, kernel_type = kernel_type)
    fit.1.p <- rd2d_lm(dat.centered[dat.centered$d == 1,], h.1, p, vce, kernel = kernel, C = C[dat.centered$d == 1],
                       varr = TRUE, kernel_type = kernel_type)
    mu.0 <- (vec %*% fit.0.p$beta)[1,1]
    mu.1 <- (vec %*% fit.1.p$beta)[1,1]

    # standardize
    if (length(h.0) == 1){
      h.0.x <- as.numeric(h.0)
      h.0.y <- as.numeric(h.0)
    }  else {
      h.0.x <- as.numeric(h.0[1])
      h.0.y <- as.numeric(h.0[2])
    }
    if (length(h.1) == 1){
      h.1.x <- as.numeric(h.1)
      h.1.y <- as.numeric(h.1)
    }  else {
      h.1.x <- as.numeric(h.1[1])
      h.1.y <- as.numeric(h.1[2])
    }

    # standard deviation
    invH.0 <- get_invH(c(h.0.x, h.0.y),p)
    se.0 <- matrix(vec, nrow = 1) %*% invH.0 %*% fit.0.p$cov.const %*% invH.0 %*% matrix(vec, ncol = 1) / (h.0.x * h.0.y)
    # se.0 <- matrix(vec, nrow = 1) %*% fit.0.p$cov.const %*% matrix(vec, ncol = 1) / (N * h.0^(2 + 2 * o))
    se.0 <- sqrt(se.0[1,1])
    invH.1 <- get_invH(c(h.1.x, h.1.y),p)
    se.1 <- matrix(vec, nrow = 1) %*% invH.1 %*% fit.1.p$cov.const %*% invH.1 %*% matrix(vec, ncol = 1) / (h.1.x * h.1.y)
    # se.1 <- matrix(vec, nrow = 1) %*% fit.1.p$cov.const %*% matrix(vec, ncol = 1) / (N * h.1^(2 + 2 * o))
    se.1 <- sqrt(se.1[1,1])

    # effective sample size
    eN.0 <- fit.0.p$eN
    eN.1 <- fit.1.p$eN

    if (ncol(hgrid.0) == 1){
      results[i,] <- c(ev[1], ev[2], h.0, h.1, mu.0, mu.1, se.0, se.1, eN.0, eN.1)
    } else {
      results[i,] <- c(ev[1], ev[2], h.0.x, h.0.y, h.1.x, h.1.y, mu.0, mu.1, se.0, se.1, eN.0, eN.1)
    }
  }

  return(results)
}

######################### Estimating Covariance matrix  ########################

rdbw2d_cov <- function(dat, eval, deriv = NULL, o = 0, p = 1, hgrid.0, hgrid.1 = NULL,
                       kernel = "epa", kernel_type = "prod", vce = "hc2", C = NULL){

  dat <- dat[,c("x.1", "x.2", "y", "d")]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]

  neval <- dim(eval)[1]
  covs <- matrix(NA, nrow = neval, ncol = neval)
  halves.0 <- list()
  halves.1 <- list()
  inds.0 <- list()
  inds.1 <- list()

  clusters <- unique(C)
  g <- length(clusters)

  # Standardize: if hgrid.0 is a vector, convert it to a m by 1 data frame
  hgrid.0 <- as.matrix(hgrid.0)
  if (!is.null(hgrid.1)) {
    hgrid.1 <- as.matrix(hgrid.1)
  }

  for (i in 1:neval){

    ev.a <- eval[i,]
    h.a.0 <- hgrid.0[i,]
    h.a.1 <- h.a.0
    if (!is.null(hgrid.1)) h.a.1 <- hgrid.1[i,]
    deriv.vec.a <- deriv[i,]

    dat.centered <- dat[, c("x.1", "x.2", "y", "d")]

    dat.centered$x.1 <- dat.centered$x.1 - ev.a$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev.a$x.2
    dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2)

    C.0 <- C[as.logical(dat.centered$d == 0)]
    C.1 <- C[as.logical(dat.centered$d == 1)]

    cov.half.consts.0 <- get_cov_half_v2(dat.centered[dat.centered$d == 0,], h.a.0, p, vce, kernel, kernel_type, C.0, clusters)
    cov.half.consts.1 <- get_cov_half_v2(dat.centered[dat.centered$d == 1,], h.a.1, p, vce, kernel, kernel_type, C.1, clusters)
    half.const.0 <- cov.half.consts.0$cov.half.const
    half.const.1 <- cov.half.consts.1$cov.half.const
    ind.0 <- cov.half.consts.0$ind
    ind.1 <- cov.half.consts.1$ind

    inds.0[[i]] <- ind.0
    inds.1[[i]] <- ind.1
    halves.0[[i]] <- half.const.0
    halves.1[[i]] <- half.const.1

  }

  for (i in 1:neval){
    for (j in i:neval){

      h.a.0 <- hgrid.0[i,]
      h.a.1 <- h.a.0
      if (!is.null(hgrid.1)) h.a.1 <- hgrid.1[i,]
      h.b.0 <- hgrid.0[j,]
      h.b.1 <- h.b.0
      if (!is.null(hgrid.1)) h.b.1 <- hgrid.1[j,]

      if (length(h.a.0) == 1){
        h.a.0.x <- h.a.0
        h.a.0.y <- h.a.0
        h.a.1.x <- h.a.1
        h.a.1.y <- h.a.1
        h.b.0.x <- h.b.0
        h.b.0.y <- h.b.0
        h.b.1.x <- h.b.1
        h.b.1.y <- h.b.1
      } else{
        h.a.0.x <- h.a.0[1]
        h.a.0.y <- h.a.0[2]
        h.a.1.x <- h.a.1[1]
        h.a.1.y <- h.a.1[2]
        h.b.0.x <- h.b.0[1]
        h.b.0.y <- h.b.0[2]
        h.b.1.x <- h.b.1[1]
        h.b.1.y <- h.b.1[2]
      }
      deriv.vec.a <- deriv[i,]
      deriv.vec.b <- deriv[j,]

      half.0.a <- halves.0[[i]]
      half.0.b <- halves.0[[j]]
      half.1.a <- halves.1[[i]]
      half.1.b <- halves.1[[j]]

      invH.a.0 <- get_invH(c(h.a.0.x,h.a.0.y),p)
      invH.a.1 <- get_invH(c(h.a.1.x,h.a.1.y),p)
      invH.b.0 <- get_invH(c(h.b.0.x,h.b.0.y),p)
      invH.b.1 <- get_invH(c(h.b.1.x,h.b.1.y),p)

      if (is.null(C)){

        ind.0.a <- inds.0[[i]]
        ind.0.b <- inds.0[[j]]
        ind.1.a <- inds.1[[i]]
        ind.1.b <- inds.1[[j]]

        ind.0 <- ind.0.a * ind.0.b
        ind.0.a <- as.logical(ind.0[ind.0.a])
        ind.0.b <- as.logical(ind.0[ind.0.b])
        ind.1 <- ind.1.a * ind.1.b
        ind.1.a <- as.logical(ind.1[ind.1.a])
        ind.1.b <- as.logical(ind.1[ind.1.b])

        cov.0 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.0 %*% t(half.0.a[ind.0.a,,drop = "FALSE"]) %*% half.0.b[ind.0.b,,drop = "FALSE"] %*% invH.b.0 %*% matrix(deriv.vec.b)
        cov.0 <- cov.0[1,1] / (sqrt(h.a.0.x * h.a.0.y * h.b.0.x * h.b.0.y))
        # cov.0 <- cov.0[1,1] / (sqrt(N * h.a^(2 + 2 * o)) * sqrt(N * h.b^(2 + 2 * o) ))
        cov.1 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.1 %*% t(half.1.a[ind.1.a,,drop = "FALSE"]) %*% half.1.b[ind.1.b,,drop = "FALSE"] %*% invH.b.1 %*% matrix(deriv.vec.b)
        cov.1 <- cov.1[1,1] / (sqrt(h.a.1.x * h.a.1.y * h.b.1.x * h.b.1.y))
        # cov.1 <- cov.1[1,1] / (sqrt(N * h.a^(2 + 2 * o)) * sqrt(N * h.b^(2 + 2 * o) ))
      }

      if (!is.null(C)){
        k <- dim(deriv)[2]
        M.0 <- matrix(0, k, k)
        M.1 <- matrix(0, k, k)
        for (l in 1:g){
          M.0 <- M.0 + matrix(half.0.a[l,], ncol = 1) %*% matrix(half.0.b[l,], nrow = 1)
          M.1 <- M.1 + matrix(half.1.a[l,], ncol = 1) %*% matrix(half.1.b[l,], nrow = 1)
        }
        cov.0 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.0 %*% M.0 %*% invH.b.0 %*% matrix(deriv.vec.b)
        cov.0 <- cov.0[1,1] /  (sqrt(h.a.0.x * h.a.0.y * h.b.0.x * h.b.0.y))
        cov.1 <- matrix(deriv.vec.a, nrow = 1) %*% invH.a.1 %*% M.1 %*% invH.b.1 %*% matrix(deriv.vec.b)
        cov.1 <- cov.1[1,1] /  (sqrt(h.a.1.x * h.a.1.y * h.b.1.x * h.b.1.y))
      }

    covs[i,j] <- cov.0 + cov.1
    covs[j,i] <- cov.0 + cov.1
    }
  }
  return(covs)
}


################################## infl ########################################

infl <- function(x, invG){
  result <- t(matrix(x, ncol = 1)) %*% invG %*% matrix(x, ncol = 1)
  return(result[1,1])
}

######################### get_cov_half (efficient) #############################

get_cov_half_v2 <- function(dat, h, p, vce = "hc1", kernel = "epa", kernel_type = "prod",
                            C = NULL, clusters = NULL){

  dat <- dat[,c("x.1", "x.2", "y", "d", "dist")]

  h <- as.vector(as.matrix(h)) # if h is data frame, convert it to a vector

  # checks
  if (kernel_type == "prod"){ # product kernel
    if (length(h) == 1){
      h <- c(h,h)
    }
  }
  else{ # radius kernel
    if (length(h) == 2){
      h <- sqrt(h[1]^2 + h[2]^2)
    }
  }

  # weights
  if (kernel_type == "prod"){
    w <- W.fun(dat$x.1/c(h[1]), kernel) * W.fun(dat$x.2/c(h[2]), kernel) / c(h[1] * h[2])
  }
  else{
    w <- W.fun(dat$dist/c(h), kernel)/c(h^2)
  }

  if (length(h) == 1){
    h.x <- h
    h.y <- h
  } else {
    h.x <- h[1]
    h.y <- h[2]
  }

  # Variance and coefficients for a linear combination of (p+1)-th derivatives

  ind <- as.logical(w > 0)

  eN <- sum(ind)

  ew <- w[ind]
  eY <- dat$y[ind]
  eC <- C[ind]

  eu <- dat[ind, c("x.1", "x.2")]
  eu$x.1 <- eu$x.1/h.x
  eu$x.2 <- eu$x.2/h.y

  eR <- as.matrix(get_basis(eu,p))

  sqrtw_R <- sqrt(ew) * eR
  sqrtw_Y <- sqrt(ew) * eY

  w_R <- ew * eR

  invG <- qrXXinv(sqrtw_R)

  invH.p <- get_invH(c(h.x,h.y),p)

  H.p <- get_H(c(h.x,h.y),p)

  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% matrix(sqrtw_Y, ncol = 1)

  resd <- abs(eY - eR %*% (H.p %*% beta))[,1]

  lambda <- function(x){infl(x, invG)}

  if (vce=="hc0") {
    w.vce <- 1
  } else if (vce=="hc1") { # TODO:Set to default
    w.vce <- sqrt(eN/(eN-factorial(p+2)/(factorial(p) * 2)))
  } else if (vce=="hc2") {
    hii <- apply(sqrtw_R, 1, lambda)
    w.vce <- sqrt(1/(1-hii))
  } else if (vce == "hc3"){
    hii <- apply(sqrtw_R, 1, lambda)
    w.vce <- 1/(1-hii)
  }

  resd <- resd * w.vce

  if (is.null(C)){

    sigma.half.const <-  (sqrt(ew) * resd * as.matrix(sqrtw_R)) * sqrt(h.x * h.y)

    cov.half.const <- sigma.half.const %*% invG
  }

  if (!is.null(C)){

    n <- length(eC)
    k <- dim(w_R)[2]

    g     <- length(clusters)
    w.w   <- ((n-1)/(n-k))*(g/(g-1))

    cov.half.const <- matrix(0,nrow = g, ncol = k)

    for (i in 1:g) {

      ind.vce <- as.logical(eC==clusters[i])
      w_R_i <- w_R[ind.vce,,drop=FALSE]
      resd_i <- resd[ind.vce]
      resd_i <- matrix(resd_i, ncol = 1)
      w_R_resd_i <- t(crossprod(w_R_i,resd_i)) * sqrt(h.x * h.y)

      cov.half.const[i,] <- as.vector(w_R_resd_i)
    }

    cov.half.const <- cov.half.const %*% invG
  }

  return(list("ind" = ind, "cov.half.const" = cov.half.const))
}

############################## Critical Value ##################################

rd2d_cval <- function(cov, rep, side="two", alpha, lp=Inf) {
  tvec <- c()

  cval <- NA
  m <- dim(cov)[1]
  cov.t <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m){
    for (j in 1:m){
      cov.t[i,j] <- cov[i,j]/sqrt(cov[i,i] * cov[j,j])
    }
  }

  sim <- mvrnorm(n = rep, mu = rep(0,m), cov.t)

  if (!is.null(side)) {
    if (side == "two") {
      if (is.infinite(lp)) {tvec <- apply(sim, c(1), function(x){max(abs(x))})}
      else                 {tvec <- apply(sim, c(1), function(x){mean(abs(x)^lp)^(1/lp)})}
    } else if (side == "left") {
      tvec <- apply(sim, c(1), max)
    } else if (side == "right") {
      tvec <- apply(sim, c(1), min)
    }
  }

  if (!is.null(side)) {
    cval <- quantile(tvec, alpha/100, na.rm=T, names = F, type=2)
  }
  return(cval)
}

############################## p Value ##################################

rd2d_pval <- function(tstat, cov, rep, side="two", lp=Inf) {
  tvec <- c()

  pval <- NA
  m <- dim(cov)[1]
  cov.t <- matrix(NA, nrow = m, ncol = m)
  for (i in 1:m){
    for (j in 1:m){
      cov.t[i,j] <- cov[i,j]/sqrt(cov[i,i] * cov[j,j])
    }
  }

  sim <- mvrnorm(n = rep, mu = rep(0,m), cov.t)

  if (!is.null(side)) {
    if (side == "two") {
      if (is.infinite(lp)) {tvec <- apply(sim, c(1), function(x){max(abs(x))})}
      else                 {tvec <- apply(sim, c(1), function(x){mean(abs(x)^lp)^(1/lp)})}
    } else if (side == "left") {
      tvec <- apply(sim, c(1), max)
    } else if (side == "right") {
      tvec <- apply(sim, c(1), min)
    }
  }

  if (!is.null(side)) {
    pval <- mean(tvec >= abs(tstat))
  }
  return(pval)
}

############################## Confidence Bands ################################

rd2d_cb <- function(mu.hat, cov.us, rep, side, alpha){

  # mu.hat: estimated (derivatives) of treatment effect
  # cov.us: estimated covariance matrix for treatment effects at all evaluation points
  # rep: number of repetitions for Gaussian simulation
  # side: "pos", "neg", or "two"
  # alpha: confidence level
  # If side == "two", returns upper and lower bounds of CI and CB
  # If side == "pos", returns upper bounds of CI and CB
  # If side == "neg", returns lower bounds of CI and CB

  se.hat <- sqrt(diag(cov.us))
  cval <- rd2d_cval(cov.us, rep = rep, side=side, alpha = alpha, lp=Inf)

  if (side == "two"){
    zval <- qnorm((alpha + 100)/ 200)
    CI.l <- mu.hat - zval * se.hat; CI.r <- mu.hat + zval * se.hat
    CB.l <- mu.hat - cval * se.hat; CB.r <- mu.hat + cval * se.hat
  }

  if (side == "left"){
    zval <- qnorm(alpha / 100)
    CI.r <- mu.hat + zval * se.hat; CI.l <- rep(-Inf, length(CI.r))
    CB.r <- mu.hat + cval * se.hat; CB.l <- rep(-Inf, length(CB.r))
  }
  if (side == "right"){
    zval <- qnorm(alpha / 100)
    CI.l <- mu.hat - zval * se.hat; CI.r <- rep(Inf, length(CI.l))
    CB.l <- mu.hat - cval * se.hat; CB.r <- rep(Inf, length(CB.l))
  }

  return(list(CI.l = CI.l, CI.r = CI.r, CB.l = CB.l, CB.r = CB.r))
}

############################# Get Basis ########################################

get_basis <- function(u,p){
  u.x.1 <- u[,1]
  u.x.2 <- u[,2]
  result <- matrix(NA, nrow = dim(u)[1], ncol = factorial(p+2)/(factorial(p) * 2))
  result[,1] <- rep(1, dim(u)[1])
  count <- 2
  if (p >= 1){
    for (j in 1:p){
      for (k in 0:j){
        result[,count] <- u.x.1^(j-k) * u.x.2^k
        count <- count + 1
      }
    }
  }
  return(result)
}

############################### Get H ##########################################

# old version with one h

# get_H <- function(h, p){
#   diags <- c()
#   for (i in 0:p){
#     diags <- c(diags, rep(h^i, i+1))
#   }
#   result <- diag(diags)
#   return(result)
# }

# new version with on h for each coordinate
get_H <- function(h,p){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }
  result <- rep(NA, factorial(p+2)/(factorial(p) * 2))
  result[1] <- 1

  if (p >= 1){
    count <- 2
    for (j in 1:p){
      for (k in 0:j){
        result[count] <- h.x.1^(j-k) * h.x.2^k
        count <- count + 1
      }
    }
  }
  result <- diag(result)
  return(result)
}

########################### Get Inverse of H ###################################

# old version with one h

# get_invH <- function(h, p){
#   diags <- c()
#   for (i in 0:p){
#     diags <- c(diags, rep(h^(-i), i+1))
#   }
#   result <- diag(diags)
#   return(result)
# }

# new version with one h for each dimension
get_invH <- function(h,p){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }
  result <- rep(NA, factorial(p+2)/(factorial(p) * 2))
  result[1] <- 1

  if (p >= 1){
    count <- 2
    for (j in 1:p){
      for (k in 0:j){
        result[count] <- 1/(h.x.1^(j-k) * h.x.2^k)
        count <- count + 1
      }
    }
  }
  result <- diag(result)
  return(result)
}

############################## lm fit ##########################################

# kernel_type = "prod" or "rad"

rd2d_lm <- function(dat, h, p, vce = "hc1", kernel = "epa", kernel_type = "prod",
                    C = NULL, varr = FALSE){

  dat <- dat[,c("x.1", "x.2", "y", "d", "dist")]

  # Variance and coefficients for a linear combination of (p+1)-th derivatives.

  h <- as.vector(as.matrix(h)) # if h is data frame, convert it to a vector

  # checks
  if (kernel_type == "prod"){ # product kernel
    if (length(h) == 1){
      h <- c(h,h)
    }
  }
  else{ # radius kernel
    if (length(h) == 2){
      h <- sqrt(h[1]^2 + h[2]^2)
    }
  }

  # weights
  if (kernel_type == "prod"){
    w <- W.fun(dat$x.1/c(h[1]), kernel) * W.fun(dat$x.2/c(h[2]), kernel) / c(h[1] * h[2])
  }
  else{
    w <- W.fun(dat$dist/c(h), kernel)/c(h^2)
  }

  if (length(h) == 1){
    h.x <- h; h.y <- h
  } else {
    h.x <- h[1]; h.y <- h[2]
  }

  ind <- as.logical(w > 0)

  eN <- sum(ind)

  ew <- w[ind]
  eY <- dat$y[ind]
  eC <- C[ind] # if C == NULL, eC == NULL.

  eu <- dat[ind, c("x.1", "x.2")]
  eu$x.1 <- eu$x.1/h.x
  eu$x.2 <- eu$x.2/h.y

  eR <- as.matrix(get_basis(eu,p))

  sqrtw_R <- sqrt(ew) * eR
  sqrtw_Y <- sqrt(ew) * eY

  w_R <- ew * eR

  invG <- qrXXinv(sqrtw_R)

  invH.p <- get_invH(c(h.x,h.y),p)

  H.p <- get_H(c(h.x,h.y),p)

  beta <- invH.p %*% invG %*% t(sqrtw_R) %*% matrix(sqrtw_Y, ncol = 1)

  cov.const <- NA

  if (varr){

    resd <- (eY - (eR %*% H.p) %*% beta)[,1]

    lambda <- function(x){infl(x, invG)}

    if (vce=="hc0") {
      w.vce = 1
    } else if (vce=="hc1") {
      w.vce = sqrt(eN/(eN-factorial(p+2)/(factorial(p) * 2)))
    } else if (vce=="hc2") {
      hii <- apply(sqrtw_R, 1, lambda)
      w.vce = sqrt(1/(1-hii))
    } else if (vce == "hc3"){
      hii <- apply(sqrtw_R, 1, lambda)
      w.vce = 1/(1-hii)
    }

    resd <- resd * w.vce

    sigma <- rd2d_vce(w_R, resd, eC, c(h.x, h.y))

    # sigma <-  t(resd * as.matrix(sqrtw_R)) %*% (ew * resd * as.matrix(sqrtw_R)) * h^2

    cov.const <- t(invG) %*% sigma %*% invG
  }

  return(list("beta" = beta, "cov.const" = cov.const, "eN" = eN))
}

###### Get coefficients of a linear combination of (p+1)-th derivatives ########

get_coeff <- function(dat.centered,vec, p,dn, kernel, kernel_type){

  dat.centered <- dat.centered[,c("x.1", "x.2", "y", "d", "dist")]
  if (kernel_type == "prod"){
    w.v <- W.fun(dat.centered$x.1/c(dn), kernel) * W.fun(dat.centered$x.2/c(dn), kernel) / c(dn * dn)
  }
  else{
    w.v <- W.fun(dat.centered$dist/c(dn), kernel)/c(dn^2)
  }

  # w.v <- W.fun(dat.centered$dist/c(dn), kernel)/c(dn^2)

  ind.v <- as.logical(w.v > 0)
  eN.v <- sum(ind.v)

  ew.v <- w.v[ind.v]
  eY.v <- dat.centered$y[ind.v]

  eu.v <- dat.centered[ind.v, c("x.1", "x.2")]
  eu.v$x.1 <- eu.v$x.1/dn
  eu.v$x.2 <- eu.v$x.2/dn

  eR.v.aug <- as.matrix(get_basis(eu.v,p+1))
  eR.v <- eR.v.aug[,1: (factorial(p+2)/(factorial(p)*2))]
  eS.v <- eR.v.aug[, (factorial(p+2)/(factorial(p)*2)+1) : (factorial(p+1+2)/(factorial(p+1)*2))]

  sqrtw_R.v <- sqrt(ew.v) * eR.v
  sqrtw_eS.v <- sqrt(ew.v) * eS.v
  sqrtw_Y.v <- sqrt(ew.v) * eY.v

  invG.v <- qrXXinv(sqrtw_R.v)

  vec.q <- matrix(vec, nrow = 1) %*% invG.v %*% t(sqrtw_R.v) %*% sqrtw_eS.v
  vec.q <- vec.q[1,]
  vec.q <- c(rep(0, factorial(p + 2)/(factorial(p)*2)), vec.q)

  return(vec.q)
}


###################### Robust Covariance Estimation ############################

# old version with one h
# rd2d_vce <- function(w_R, resd, eC, h){
#   n <- length(eC)
#   k <- dim(w_R)[2]
#   M <- matrix(0, nrow = k, ncol = k)
#   if (is.null(eC)){
#     w.w <- 1
#     M <-  crossprod(resd * as.matrix(w_R)) * h^2
#   }
#   else{
#     clusters = unique(eC)
#     g     = length(clusters)
#     w.w =((n-1)/(n-k))*(g/(g-1))
#     for (i in 1:g) {
#       ind=eC==clusters[i]
#       w_R_i = w_R[ind,,drop=FALSE]
#       resd_i = resd[ind]
#       resd_i <- matrix(resd_i, ncol = 1)
#       w_R_resd_i = t(crossprod(w_R_i,resd_i))
#       M = M + crossprod(w_R_resd_i,w_R_resd_i) * h^2
#     }
#   }
#   return(M * w.w)
# }

# new version with one h for each coordinate
rd2d_vce <- function(w_R, resd, eC, h){
  if (length(h) == 1){
    h.x.1 <- h
    h.x.2 <- h
  } else {
    h.x.1 <- h[1]
    h.x.2 <- h[2]
  }

  n <- length(eC)
  k <- dim(w_R)[2]
  M <- matrix(0, nrow = k, ncol = k)
  if (is.null(eC)){
    w.w <- 1
    M <-  crossprod(resd * as.matrix(w_R)) * h.x.1 * h.x.2
  }
  else{
    clusters = unique(eC)
    g     = length(clusters)
    w.w =((n-1)/(n-k))*(g/(g-1))
    for (i in 1:g) {
      ind=eC==clusters[i]
      # w_R_i = w_R[ind,,drop=FALSE]
      # Attempt to subset w_R with ind
      w_R_i <- tryCatch(
        {
          w_R[ind, , drop = FALSE]
        },
        error = function(e) {
          cat("Error: ", conditionMessage(e), "\n")
          cat("Dimensions of w_R: ", paste(dim(w_R), collapse = " x "), "\n")
          cat("Length of ind: ", length(ind), "\n")
          stop("Exiting due to the above error.")
        }
      )
      resd_i = resd[ind]
      resd_i <- matrix(resd_i, ncol = 1)
      w_R_resd_i = t(crossprod(w_R_i,resd_i))
      M = M + crossprod(w_R_resd_i,w_R_resd_i) * h.x.1 * h.x.2
    }
  }
  return(M * w.w)
}

############################## rd2d unique #####################################

rd2d_unique <- function(dat){

  dat <- dat[,c("x.1", "x.2", "y", "d")]

  ord <- order(dat[,1], dat[,2])

  dat <- dat[ord,]

  N <- dim(dat)[1]

  # if x has one or no element
  if (N == 0) return(list(unique = NULL, freq = c(), index = c()))
  if (N == 1) return(list(unique = dat, freq = 1, index = 1))

  # else
  uniqueIndex <- c(c(dat[2:N, 1] != dat[1:(N-1),1] | dat[2:N, 2] != dat[1:(N-1),2]), TRUE)
  unique <- dat[uniqueIndex,]
  nUnique <- dim(unique)[1]

  # all are distinct
  if (nUnique == N) return(list(unique=unique, freq=rep(1,N), index=1:N))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=N, index=N))

  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1

  return(list(unique=unique, freq=freq, index=(1:N)[uniqueIndex]))
}

############################# output formatting ################################

# TODO: merge kernel and kernel.type

# TODO: merge cluster and vce
