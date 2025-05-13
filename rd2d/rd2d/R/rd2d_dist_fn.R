# library(MASS)
# # library(tidyr)
# # library("lmtest")
# # library("sandwich")
# library(expm)
#
# source("rd2d_fn_v2.R")

############################ Rule of Thumb #####################################

# Use Scott's Rule: n^{-1/(d+4)} * sample variance

rdbw2d_dist_rot <- function(dist, kernel.type){

  mu2K.squared <- NA
  l2K.squared <- NA

  if (kernel.type == "Epanechnikov"){mu2K.squared <- 1/6; l2K.squared <- 4/(3 * pi)}
  if (kernel.type == "Triangular"){mu2K.squared <- 3/20; l2K.squared <- 3/(2 * pi)}
  if (kernel.type == "Uniform"){mu2K.squared <- 1/4; l2K.squared <- 1/(pi)}
  if (kernel.type == "Gaussian"){mu2K.squared <- 1; l2K.squared <- 1/(4 * pi)}

  # Estimate sample variance.

  N <- length(dist)

  var.hat <- 1 / 2 * mean(dist^2)
  cov.matrix <- diag(c(var.hat, var.hat))

  D <- 2

  trace.const <- 1/( 2^(D+2) * pi^(D/2) * det(sqrtm(cov.matrix)) ) * ( 2*sum(diag(ginv(cov.matrix, 1e-20) %*% ginv(cov.matrix, 1e-20)))
                                                                       + (sum(diag(ginv(cov.matrix, 1e-20))))^2 )

  hROT <-( (D * l2K.squared) / (N * mu2K.squared * trace.const) )^(1/(4+D))

  return(hROT)
}

############# One Step Bandwidth Selection for Distance Based Estimator ################

rdbw2d_dist_bw <- function(Y, D, p = 1, kernel, deriv.vec = NULL,
                  rot = NULL, vce = "hc0", C = NULL, bwcheck = 50 + p + 1,
                  scaleregul = 1, cqt = 0.5, verbose = FALSE){

  # Data Cleaning

  neval <- ncol(D)
  N <- length(Y)
  if (!is.null(C)){
    M <- length(unique(C))
  } else {
    M <- N
  }

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  # Loop over points of evaluations

  results <- data.frame(matrix(NA, ncol = 14, nrow = neval))
  colnames(results) <- c('h.0', 'h.1', 'b.0', 'b.1', 'v.0', 'v.1','r.0','r.1','eN.0', 'eN.1', 'bw.min.0', 'bw.min.1', 'bw.max.0', 'bw.max.1')

  for (i in 1:neval){

    dist <- D[,i]
    d <- (dist >= 0)
    dist <- abs(dist)

    if (!is.null(rot)){
      dn <- rot[i]
    } else {
      dn <- rdbw2d_dist_rot(dist, kernel.type)
    }

    vec <- deriv.vec[i,]

    # Estimate (p+1)-th derivatives using 50% of data

    dat.nearby <- cbind(Y,d,dist)
    dat.nearby <- as.data.frame(dat.nearby)
    colnames(dat.nearby) <- c("y", "d", "dist")

    thrshd.0 <- quantile((dat.nearby[dat.nearby$d == 0,])$dist, cqt)
    thrshd.1 <- quantile((dat.nearby[dat.nearby$d == 1,])$dist, cqt)

    dat.nearby.nm <- c()
    j <- 2
    count <- 1

    while (j <= p+1){
      dat.nearby <- cbind(dat.nearby, (dat.nearby$dist)^j )
      dat.nearby.nm <- c(dat.nearby.nm, paste("V", count, sep = ''))
      count <- count + 1
      j <- j+1
    }
    colnames(dat.nearby) <- c(c("y", "d", "dist"),dat.nearby.nm)

    count <- count - 1
    fml <-"y ~ dist"
    if (count >= 1){
      for (k in 1:count){
        fml <- paste(fml," + ","V",k,sep = "")
      }
    }
    fml <- as.formula(fml)

    # Fit lp for (p + 1)th derivatives

    model.deriv.0 <- lm(fml, data = dat.nearby[(dat.nearby$d == FALSE) & (dat.nearby$dist <= thrshd.0),])
    model.deriv.1 <- lm(fml, data = dat.nearby[(dat.nearby$d == TRUE) & (dat.nearby$dist <= thrshd.1),])

    deriv.ppls1.0 <- model.deriv.0$coefficients[length(model.deriv.0$coefficients)]
    deriv.ppls1.1 <- model.deriv.1$coefficients[length(model.deriv.1$coefficients)]

    if (verbose){
      print("deriv.ppls1.0"); print(deriv.ppls1.0)
      print("deriv.ppls1.1"); print(deriv.ppls1.1)
    }

    # Take effective sample
    dn.0 <- dn; dn.1 <- dn

    bw.min.0 <- NA
    bw.min.1 <- NA
    bw.max.0 <- NA
    bw.max.1 <- NA

    if (!is.null(bwcheck)) { # Bandwidth restrictions
      sorted.0 <- sort(dat.nearby[dat.nearby$d == FALSE,]$dist)
      sorted.1 <- sort(dat.nearby[dat.nearby$d == TRUE,]$dist)
      bw.min.0   <- sorted.0[min(bwcheck,length(sorted.0))]
      bw.min.1   <- sorted.1[min(bwcheck,length(sorted.1))]
      bw.max.0   <- sorted.0[length(sorted.0)]
      bw.max.1   <- sorted.1[length(sorted.1)]
      dn.0     <- max(dn, bw.min.0)
      dn.1     <- max(dn, bw.min.1)
      dn.0     <- min(dn, bw.max.0)
      dn.1     <- min(dn, bw.max.1)
    }

    w.0   <- W.fun(dat.nearby[dat.nearby$d == FALSE,]$dist/dn.0, kernel)/(dn.0^2)
    w.1   <- W.fun(dat.nearby[dat.nearby$d == TRUE,]$dist/dn.1, kernel)/(dn.1^2)
    ind.0 <- (w.0 > 0); ind.1 <- (w.1 > 0)
    eN.0 <- sum(ind.0); eN.1 <- sum(ind.1)

    eX.0  <- dat.nearby[dat.nearby$d==FALSE,][c("dist")][ind.0,] # Distance instead of coordinates
    eX.1  <- dat.nearby[dat.nearby$d==TRUE,][c("dist")][ind.1,]
    eY.0 <- dat.nearby[dat.nearby$d==FALSE,][c("y")][ind.0,] # df["A"] will always return a data frame
    eY.1 <- dat.nearby[dat.nearby$d==TRUE,][c("y")][ind.1,]
    eW.0 <- w.0[ind.0]
    eW.1 <- w.1[ind.1]
    u.0   <- eX.0
    u.1   <- eX.1

    # -- Sd from residuals of local polynomial fit

    sd.0 <- eY.0 - predict(model.deriv.0, newdata = dat.nearby[dat.nearby$d==FALSE,][ind.0,])
    sd.1 <- eY.1 - predict(model.deriv.1, newdata = dat.nearby[dat.nearby$d==TRUE,][ind.1,])

    # Compute bread, meat and half-bread matrices for p-th degree model

    R.0.p <- as.matrix(get_basis_v1(u.0/dn.0,p))
    R.1.p <- as.matrix(get_basis_v1(u.1/dn.1,p))
    inv.gamma0.p <- ginv(crossprod(sqrt(unname(eW.0)) * as.matrix(R.0.p)), 1e-20) # Bread matrix
    inv.gamma1.p <- ginv(crossprod(sqrt(unname(eW.1)) * as.matrix(R.1.p)), 1e-20) # Bread matrix
    sigma0.half.p <- eW.0 * as.matrix(R.0.p)
    sigma1.half.p <- eW.1 * as.matrix(R.1.p) # Half-bread matrix

    # vce methods

    lambda.0 <- function(x){infl(x, inv.gamma0.p)}
    lambda.1 <- function(x){infl(x, inv.gamma1.p)}

    sqrtw_R.0 <- sqrt(eW.0) * as.matrix(R.0.p)
    sqrtw_R.1 <- sqrt(eW.1) * as.matrix(R.1.p)

    if (vce=="hc0") {
      w.vce.0 = 1
      w.vce.1 = 1
    } else if (vce=="hc1") {
      w.vce.0 = sqrt(eN.0/(eN.0-p-1))
      w.vce.1 = sqrt(eN.1/(eN.1-p-1))
    } else if (vce=="hc2") {
      hii.0 <- apply(sqrtw_R.0, 1, lambda.0)
      hii.1 <- apply(sqrtw_R.1, 1, lambda.1)
      w.vce.0 = sqrt(1/(1-hii.0))
      w.vce.1 = sqrt(1/(1-hii.1))
    } else if (vce=="hc3"){
      hii.0 <- apply(sqrtw_R.0, 1, lambda.0)
      hii.1 <- apply(sqrtw_R.1, 1, lambda.1)
      w.vce.0 = 1/(1-hii.0)
      w.vce.1 = 1/(1-hii.1)
    }
    sd.0 <- sd.0 * w.vce.0
    sd.1 <- sd.1 * w.vce.1

    # (cluster-robust) meat matrix estimation

    eC.0 <- C[dat.nearby$d == FALSE][ind.0]
    eC.1 <- C[dat.nearby$d == TRUE][ind.1]
    sigma.0.p <- rd2d_vce(sigma0.half.p, sd.0, eC.0, dn.0)
    sigma.1.p <- rd2d_vce(sigma1.half.p, sd.1, eC.1, dn.1)

    # sigma.0.p <-  t(sd.0 * sigma0.half.p) %*% (sd.0 * sigma0.half.p) * eN.0 / (eN.0 - (p+1)) * dn^2 # Meat matrices
    # sigma.1.p <-  t(sd.1 * sigma1.half.p) %*% (sd.1 * sigma1.half.p) * eN.1 / (eN.1 - (p+1)) * dn^2 # Meat matrices

    # Compute the coefficients for linear combination of (p+1)-th derivatives

    pmatrix.0 <- matrix(NA,nrow = eN.0, ncol = 1)
    pmatrix.1 <- matrix(NA, nrow = eN.1, ncol = 1)
    pmatrix.0[,1] <- (eX.0/dn)^(p+1) * eW.0
    pmatrix.1[,1] <- (eX.1/dn)^(p+1) * eW.1
    pmatrix.0 <- t(R.0.p) %*% pmatrix.0
    pmatrix.1 <- t(R.1.p) %*% pmatrix.1
    dmm <- p+2
    coeff.0 <- rep(0, dmm)
    coeff.1 <- rep(0, dmm)
    coeff.0[dmm] <- as.vector( t(as.matrix(vec)) %*% inv.gamma0.p %*% pmatrix.0 )
    coeff.1[dmm] <- as.vector( t(as.matrix(vec)) %*% inv.gamma1.p %*% pmatrix.1 )

    # Compute bias constants for lp using p-th order model

    B.p.0 <- coeff.0 %*% model.deriv.0$coefficients# -- coefficient using dn
    B.p.1 <- coeff.1 %*% model.deriv.1$coefficients# -- coefficient using dn

    # Compute regularization terms for lp using (p+1)-th order model

    R.0 <- coeff.0[dmm]^2 * summary(model.deriv.0)$coefficients[dmm,2]^2
    R.1 <- coeff.1[dmm]^2 * summary(model.deriv.1)$coefficients[dmm,2]^2

    if (verbose){
      print("model.deriv.0$coefficients = ")
      print(model.deriv.0$coefficients)
      print("model.deriv.1$coefficients = ")
      print(model.deriv.1$coefficients)
    }

    # Compute variance constants for lp using p-th order model

    # -- inv.gamma0 and sigma0.half using dn.0

    V.p.0 <- t(as.matrix(vec)) %*% inv.gamma0.p %*% sigma.0.p %*% inv.gamma0.p %*% as.matrix(vec)
    V.p.1 <- t(as.matrix(vec)) %*% inv.gamma1.p %*% sigma.1.p %*% inv.gamma1.p %*% as.matrix(vec)

    # Optimal bandwidth for treatment effect using p-th order model

    hn <- (2 * (V.p.0  + V.p.1) / ( (2 * p + 2) * ((B.p.0 -  B.p.1)^2 + scaleregul * R.0 + scaleregul * R.1)) )^(1/(2 * p + 4))

    if (verbose){
      print(paste("V.p.0 = ", V.p.0, ", V.p.1 = ", V.p.1, sep = ""))
      print("coeff.0 = "); print(coeff.0)
      print("coeff.1 = "); print(coeff.1)
    }
    results[i,] <- c(hn, hn, B.p.0, B.p.1, V.p.0, V.p.1, R.0,R.1,eN.0, eN.1, bw.min.0, bw.min.1, bw.max.0, bw.max.1)
  }

  return(results)
}

############################### rd2d_dist_fit ##################################

rd2d_dist_fit <- function(Y, D, h, p, b, kernel, vce, bwcheck, masspoints, C, cbands = TRUE){

  neval <- ncol(D)

  if (!is.null(b)){
    eval <- as.data.frame(b)
  } else {
    eval <- matrix(NA, nrow = neval, ncol = 2)
    eval <- as.data.frame(eval)
  }
  colnames(eval) <- c("x.1", "x.2")

  hgrid <- h[,1]
  hgrid.1 <- h[,2]
  neval <- ncol(D)
  d <- (D[,1] >= 0)
  N <- length(Y)
  N.1 <- sum(d)
  N.0 <- N - N.1
  Estimate=matrix(NA,neval,10)
  Estimate = as.data.frame(Estimate)
  colnames(Estimate)=c("b1", "b2", "h0", "h1", "N0","N1", "mu0","mu1","se0","se1")

  # to store matrices

  Indicators.0 <- list()
  Indicators.1 <- list()
  inv.designs.0 <- list()
  inv.designs.1 <- list()
  sig.halfs.0 <- list()
  sig.halfs.1 <- list()
  resd.0 <- list()
  resd.1 <- list()

  # Check for mass points

  M.vec <- rep(N, neval)
  M.0.vec <- rep(N.0, neval)
  M.1.vec <- rep(N.1, neval)
  is_mass_point <- 0
  if (masspoints == "check" | masspoints == "adjust"){
    for (j in 1:ncol(D)){
      dist <- D[,j]
      unique.const <- rd2d_dist_unique(dist)
      unique <- unique.const$unique
      M.0 <- sum(unique <= 0)
      M.1 <- sum(unique > 0)
      M <- M.0 + M.1
      mass <- 1 - M / N
      M.vec[j] <- M
      M.0.vec[j] <- M.0
      M.1.vec[j] <- M.1
      if (mass >= 0.2){is_mass_point <- 1}
    }
    if (is_mass_point > 0){
      warning("Mass points detected in the running variables.")
      if (masspoints == "check") warning("Try using option masspoints=adjust.")
      if (is.null(bwcheck) & (masspoints == "check" | masspoints == "adjust")) bwcheck <- 50 + p + 1
    }
  }

  for (i in 1:neval) {

    b1 <- eval[i,1]
    b2 <- eval[i,2]
    dist <- D[,i]
    d <- (dist >= 0)
    dist <- abs(dist)

    dat.nearby <- cbind(Y,d,dist)
    dat.nearby <- as.data.frame(dat.nearby)
    colnames(dat.nearby) <- c("y", "d", "dist")

    h.0 <- hgrid[i]
    h.1 <- hgrid.1[i]

    M.0 <- M.0.vec[i]
    M.1 <- M.1.vec[i]

    # bandwidth restrictions
    bw.min.0 <- NA
    bw.min.1 <- NA
    bw.max.0 <- NA
    bw.max.1 <- NA

    if (!is.null(bwcheck)) { # Bandwidth restrictions
      sorted.0 <- sort(dat.nearby[dat.nearby$d == FALSE,]$dist)
      sorted.1 <- sort(dat.nearby[dat.nearby$d == TRUE,]$dist)
      bw.min.0   <- sorted.0[min(bwcheck,length(sorted.0))]
      bw.min.1   <- sorted.1[min(bwcheck,length(sorted.1))]
      bw.max.0   <- sorted.0[length(sorted.0)]
      bw.max.1   <- sorted.1[length(sorted.1)]
      h.0     <- max(h.0, bw.min.0)
      h.1     <- max(h.1, bw.min.1)
      h.0     <- min(h.0, bw.max.0)
      h.1     <- min(h.1, bw.max.1)
    }

    temp.data.control <- dat.nearby[dat.nearby$d==0,]
    temp.data.control$w <- W.fun(temp.data.control$dist/h.0, kernel = kernel)/c(h.0^2)
    temp.data.treated <- dat.nearby[dat.nearby$d==1,]
    temp.data.treated$w <- W.fun(temp.data.treated$dist/h.1, kernel = kernel)/c(h.1^2)

    w.0   <- temp.data.control$w
    w.1   <- temp.data.treated$w
    ind.0 <- (w.0 > 0)
    ind.1 <- (w.1 > 0)

    eN.0  <- sum(ind.0); eN.1 <- sum(ind.1)
    eY.0  <- temp.data.control$y[ind.0]; eY.1  <- temp.data.treated$y[ind.1]
    u.0  <- temp.data.control$dist[ind.0]; u.1  <- temp.data.treated$dist[ind.1]
    eW.0 <- w.0[ind.0]; eW.1 <- w.1[ind.1]
    eC.0 <- C[d == FALSE][ind.0]; eC.1 <- C[d == TRUE][ind.1]

    R.0 <- matrix(NA, nrow = eN.0, ncol = p+1); R.1 <- matrix(NA, nrow = eN.1, ncol = p+1)
    for (j in 1:(p+1)){
      R.0[,j] <- (u.0 / h.0)^(j-1); R.1[,j] <- (u.1 /h.1)^(j-1)
    }

    inv.gamma0 <- ginv(crossprod(sqrt(unname(eW.0)) * as.matrix(R.0)), 1e-20)
    inv.gamma1 <- ginv(crossprod(sqrt(unname(eW.1)) * as.matrix(R.1)), 1e-20)

    sigma0.half <- eW.0 * as.matrix(R.0); sigma1.half <- eW.1 * as.matrix(R.1)

    XtWY.0 <- t(matrix(eY.0 * unname(eW.0), nrow = 1) %*% as.matrix(R.0)); XtWY.1 <-t(matrix(eY.1 * unname(eW.1), nrow = 1) %*% as.matrix(R.1))
    hbeta.0 <- inv.gamma0 %*% XtWY.0; hbeta.1 <- inv.gamma1 %*% XtWY.1
    mu0 <- hbeta.0[1];  mu1 <- hbeta.1[1]
    invH.0 <- get_invH_dist(h.0,p); invH.1 <- get_invH_dist(h.1,p)
    res0 <- abs(eY.0 - as.matrix(R.0) %*% (invH.0 %*% hbeta.0))[,1]
    res1 <- abs(eY.1 - as.matrix(R.1) %*% (invH.1 %*% hbeta.1))[,1]

    sqrtw_R.0 <- sqrt(eW.0) * as.matrix(R.0)
    sqrtw_R.1 <- sqrt(eW.1) * as.matrix(R.1)

    ############################ start: get_cov_half ###########################
    if (is.null(C)){
      cov.half.const.0 <- res0 * sigma0.half %*% inv.gamma0
      cov.half.const.1 <- res1 * sigma1.half %*% inv.gamma1
    } else {
      k <- p + 1
      clusters <- unique(C)
      g <- length(clusters)
      n.0 <- length(eC.0); n.1 <- length(eC.1)
      w.w.0   <- ((n.0-1)/(n.0-k))*(g/(g-1))
      w.w.1   <- ((n.1-1)/(n.1-k))*(g/(g-1))
      cov.half.const.0 <- matrix(0,nrow = g, ncol = k)
      cov.half.const.1 <- matrix(0,nrow = g, ncol = k)

      for (l in 1:g) {
        ind.vce.0 <- as.logical(eC.0==clusters[l]); ind.vce.1 <- as.logical(eC.1==clusters[l])
        w_R_i.0 <- sigma0.half[ind.vce.0,,drop=FALSE]; w_R_i.1 <- sigma1.half[ind.vce.1,,drop=FALSE]
        resd_i.0 <- res0[ind.vce.0]; resd_i.1 <- res1[ind.vce.1]
        resd_i.0 <- matrix(resd_i.0, ncol = 1); resd_i.1 <- matrix(resd_i.1, ncol = 1)
        w_R_resd_i.0 <- t(crossprod(w_R_i.0,resd_i.0)); w_R_resd_i.1 <- t(crossprod(w_R_i.1,resd_i.1))
        cov.half.const.0[l,] <- as.vector(w_R_resd_i.0); cov.half.const.1[l,] <- as.vector(w_R_resd_i.1)
      }
      cov.half.const.0 <- cov.half.const.0 %*% inv.gamma0
      cov.half.const.1 <- cov.half.const.1 %*% inv.gamma1
    }
    ############################ end: get_cov_half #############################

    # standard deviation
    lambda.0 <- function(x){infl(x, inv.gamma0)}
    lambda.1 <- function(x){infl(x, inv.gamma1)}

    if (vce=="hc0") {
      w.vce.0 = 1
      w.vce.1 = 1
    } else if (vce=="hc1") {
      w.vce.0 = sqrt(eN.0/(eN.0- p - 1))
      w.vce.1 = sqrt(eN.1/(eN.1- p - 1))
    } else if (vce=="hc2") {
      hii.0 <- apply(sqrtw_R.0, 1, lambda.0)
      w.vce.0 = sqrt(1/(1-hii.0))
      hii.1 <- apply(sqrtw_R.1, 1, lambda.1)
      w.vce.1 = sqrt(1/(1-hii.1))
    } else if (vce == "hc3"){
      hii.0 <- apply(sqrtw_R.0, 1, lambda.0)
      w.vce.0 = 1/(1-hii.0)
      hii.1 <- apply(sqrtw_R.1, 1, lambda.1)
      w.vce.1 = 1/(1-hii.1)
    }

    res0 <- res0 * w.vce.0
    res1 <- res1 * w.vce.1

    sigma.0 <- rd2d_vce(sigma0.half, res0, eC.0, c(h.0, h.0))
    sigma.1 <- rd2d_vce(sigma1.half, res1, eC.1, c(h.1, h.1))

    cov.const.0 <- t(inv.gamma0) %*% sigma.0 %*% inv.gamma0
    cov.const.1 <- t(inv.gamma1) %*% sigma.1 %*% inv.gamma1

    vec <- rep(0,p+1)
    vec[1] <- 1
    se0 <- matrix(vec, nrow = 1) %*% cov.const.0 %*% matrix(vec, ncol = 1) / (h.0^2)
    se0 <- sqrt(se0[1,1])
    se1 <- matrix(vec, nrow = 1) %*% cov.const.1 %*% matrix(vec, ncol = 1) / (h.1^2)
    se1 <- sqrt(se1[1,1])

    if (cbands){
      Indicators.0[[i]] <- ind.0
      Indicators.1[[i]] <- ind.1
      inv.designs.0[[i]] <- inv.gamma0
      inv.designs.1[[i]] <- inv.gamma1
      sig.halfs.0[[i]] <- cov.half.const.0
      sig.halfs.1[[i]] <- cov.half.const.1
      resd.0[[i]] <- res0
      resd.1[[i]] <- res1
    }

    Estimate[i,] <- c(b1, b2, h.0, h.1, eN.0, eN.1, mu0, mu1, se0, se1)
  }

  return(list("Estimate" = Estimate, "Indicators.0" = Indicators.0, "Indicators.1" = Indicators.1,
              "inv.designs.0" = inv.designs.0, "inv.designs.1" = inv.designs.1, "sig.halfs.0" = sig.halfs.0,
              "sig.halfs.1" = sig.halfs.1, "resd.0" = resd.0, "resd.1" = resd.1,
              "M.vec" = M.vec, "M.0.vec" = M.0.vec, "M.1.vec" = M.1.vec))
}

########################### Standardization ####################################

# Standardize input data and point of evaluations.

lpgeo_bwselect_std <- function(x, eval){

  scale.1 <- sd(x$x.1)
  scale.2 <- sd(x$x.2)

  x.std <- x
  eval.std <- eval
  x.std$x.1 <- x.std$x.1 / scale.1
  x.std$x.2 <- x.std$x.2 / scale.2
  eval.std$x.1 <- eval.std$x.1 / scale.1
  eval.std$x.2 <- eval.std$x.2 / scale.2

  return(list(x.std = x.std, eval.std = eval.std, scale.1 = scale.1, scale.2 = scale.2))
}

###################################### Misc ####################################

pow <- function(vec,p){
  # if input is a vecotr
  x <- vec[1]
  y <- vec[2]
  # vec: vector
  # p: power to raise
  # d: dimension of vec
  res <- c(1)
  if (p>= 1){
    for (i in 1:p){
      for (j in 0:i){
        res <- c(res, x^(i-j) * y^(j))
      }
    }
  }
  return(res)
}

# pow.d
pow.d <- function(x,p){
  # if input is a numeric
  res <- rep(NA,p+1)
  res[1] <- 1
  if (p >= 1){
    for (i in 1:p){
      res[i+1] <- x^i
    }
  }
  return(res)
}

get_basis_v1 <- function(u,p) {
  u <- as.matrix(u)
  if (dim(u)[2] == 1){
    pow.temp <- function(v){ return(pow.d(v,p))}
  } else {
    pow.temp <- function(v){ return(pow(v,p)) }
  }
  tmp <- apply(u,c(1),pow.temp,simplify = FALSE)
  res <- data.frame(do.call(rbind,as.matrix(tmp)))
  return(res)
}

rd2d_dist_unique <- function(dist){

  ord <- order(dist)

  dist <- dist[ord]

  N <- length(dist)

  # if x has one or no element
  if (N == 0) return(list(unique = NULL, freq = c(), index = c()))
  if (N == 1) return(list(unique = dist, freq = 1, index = 1))

  # else
  uniqueIndex <- c(c(dist[2:N] != dist[1:(N-1)]), TRUE)
  unique <- dist[uniqueIndex]
  nUnique <- length(unique)

  # all are distinct
  if (nUnique == N) return(list(unique=unique, freq=rep(1,N), index=1:N))
  # all are the same
  if (nUnique == 1) return(list(unique=unique, freq=N, index=N))

  # otherwise
  freq <- (cumsum(!uniqueIndex))[uniqueIndex]
  freq <- freq - c(0, freq[1:(nUnique-1)]) + 1

  return(list(unique=unique, freq=freq, index=(1:N)[uniqueIndex]))
}

get_invH_dist <- function(h,p){
  result <- rep(NA, p+1)
  result[1] <- 1

  if (p >= 1){
    for (j in 1:p){
      result[j+1] <- 1/h^j
    }
  }
  result <- diag(result)
  return(result)
}


