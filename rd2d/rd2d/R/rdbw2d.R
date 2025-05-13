# This file is for bandwidth selection for local polynomial estimator based on 2d location data.

#' @title Bandwidth Selection for 2D Local Polynomial RD Design
#'
#' @description
#' \code{rdbw2d} implements bandwidth selector for bivariate local polynomial boundary regression discontinuity (RD) point estimators with robust bias-corrected pointwise confidence intervals and 
#' uniform confidence bands, developed in Cattaneo, Titiunik and Yu (2025). 
#' 
#' Companion commands are: \code{rd2d} for point estimation and inference procedures.
#' 
#' For other packages of RD designs, visit
#' <https://rdpackages.github.io/>
#'
#' @param Y Dependent variable; a numeric vector of length \eqn{N}, where \eqn{N} is the sample size.
#' @param X Bivariate running variable (a.k.a score variable); a numeric matrix or data frame of dimension \eqn{N \times 2}, with each row \eqn{\mathbf{X}_i = (X_{1i}, X_{2i})}.
#' @param t Treatment indicator; a logical or binary vector indicating treatment assignment (\eqn{t_i = 1} if treated, \eqn{t_i = 0} otherwise).
#' @param b Evaluation points; a matrix or data frame specifying boundary points \eqn{\mathbf{b}_j = (b_{1j}, b_{2j})}, of dimension \eqn{J \times 2}.
#' @param p Polynomial order of local polynomial estimator.
#' @param deriv The order of the derivatives of the regression functions to be estimated; a numeric vector of length 2 specifying the number of derivatives in each coordinate (e.g., \eqn{c(1,2)} corresponds to \eqn{\partial_1 \partial_2^2}).
#' @param tangvec Tangent vectors; a matrix or data frame of dimension \eqn{J \times 2} specifying directional derivatives. Overrides \code{deriv} if provided.
#' @param kernel Kernel function to use. Options are \code{"unif"}, \code{"uniform"} (uniform), \code{"triag"}, \code{"triangular"} (triangular, default), and \code{"epan"}, \code{"epanechnikov"} (Epanechnikov).
#' @param kernel_type Kernel structure. Either \code{"prod"} for product kernels or \code{"rad"} for radial kernels.
#' @param bwselect Bandwidth selection strategy. Options: 
#' \itemize{
#' \item \code{"mserd"}. One common MSE-optimal bandwidth selector for the boundary RD treatment effect estimator for each evaluation point (default). 
#' \item \code{"imserd"}. IMSE-optimal bandwidth selector for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"msetwo"}. Two different MSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator for each evaluation point. 
#' \item \code{"imsetwo"}. Two IMSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"user provided"}. User-provided bandwidths. If \code{h} is not \code{NULL}, then \code{bwselect} is overwritten to \code{"user provided"}.
#' }
#' @param method Bandwidth selection method for bias estimator based on local polynomials. Either \code{"dpi"} (default) for data-driven plug-in MSE optimal bandwidth selector or \code{"rot"} for rule-of-thumb bandwidth selector.
#' @param vce Variance-covariance estimation method. Options are:
#' \itemize{
#' \item \code{"hc0"}: heteroskedasticity-robust plug-in residual variance estimator without small-sample adjustment.
#' \item \code{"hc1"}: heteroskedasticity-robust plug-in residual variance estimator with HC1 small-sample adjustment (default).
#' \item \code{"hc2"}: heteroskedasticity-robust plug-in residual variance estimator with HC2 adjustment.
#' \item \code{"hc3"}: heteroskedasticity-robust plug-in residual variance estimator with HC3 adjustment.
#' }
#' Default is \code{"hc1"}.
#' @param bwcheck If a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least \code{bwcheck} unique observations are used. Default is \code{50 + p + 1}.
#' @param masspoints Handling of mass points in the running variable. Options are:
#' \itemize{
#' \item \code{"check"}: detects presence of mass points and reports the number of unique observations (default).
#' \item \code{"adjust"}: adjusts preliminary bandwidths to ensure a minimum number of unique observations within each side of the cutoff.
#' \item \code{"off"}: ignores presence of mass points.
#' }
#' @param C Cluster ID variable used for cluster-robust variance estimation with degrees-of-freedom weights.Default is \code{C = NULL}.
#' @param scaleregul Scaling factor for the regularization term in bandwidth selection. Default is 3.
#' @param scalebiascrct Scaling factor used for bias correction based on higher order expansions. Default is 1.
#' @param stdvars Logical. If TRUE, the running variables \eqn{X_{1i}} and \eqn{X_{2i}} are standardized before computing the bandwidths. Default is \code{FALSE}. Standardization only affects automatic bandwidth selection if bandwidths are not manually provided via \code{h}.
#'
#' @return A list of class \code{"rdbw2d"} containing:
#' \describe{
#'   \item{\code{bws}}{Data frame of estimated bandwidths for each evaluation point:
#'     \describe{
#'       \item{\code{b1}}{First coordinate of the evaluation point.}
#'       \item{\code{b2}}{Second coordinate of the evaluation point.}
#'       \item{\code{h01}}{Estimated bandwidth for \eqn{X_{1i}} in the control group (\eqn{\mathcal{A}_0}).}
#'       \item{\code{h02}}{Estimated bandwidth for \eqn{X_{2i}} in the control group (\eqn{\mathcal{A}_0}).}
#'       \item{\code{h11}}{Estimated bandwidth for \eqn{X_{1i}} in the treatment group (\eqn{\mathcal{A}_1}).}
#'       \item{\code{h12}}{Estimated bandwidth for \eqn{X_{2i}} in the treatment group (\eqn{\mathcal{A}_1}).}
#'     }
#'   }
#'   \item{\code{mseconsts}}{Data frame of intermediate quantities used in bandwidth calculation:
#'     \describe{
#'       \item{\code{Nh0}}{Effective sample size for the control group \eqn{\mathcal{A}_0}.}
#'       \item{\code{Nh1}}{Effective sample size for the treatment group \eqn{\mathcal{A}_1}.}
#'       \item{\code{bias.0}}{Bias constant estimate for the control group.}
#'       \item{\code{bias.1}}{Bias constant estimate for the treatment group.}
#'       \item{\code{var.0}}{Variance constant estimate for the control group.}
#'       \item{\code{var.1}}{Variance constant estimate for the treatment group.}
#'       \item{\code{reg.bias.0}}{Bias correction adjustment for the control group.}
#'       \item{\code{reg.bias.1}}{Bias correction adjustment for the treatment group.}
#'       \item{\code{reg.var.0}}{Variance of the bias estimate for the control group.}
#'       \item{\code{reg.var.1}}{Variance of the bias estimate for the treatment group.}
#'     }
#'   }
#'   \item{\code{opt}}{List containing:
#'     \describe{
#'       \item{\code{p}}{Polynomial order used for estimation.}
#'       \item{\code{kernel}}{Kernel function used.}
#'       \item{\code{kernel_type}}{Type of kernel (product or radial).}
#'       \item{\code{stdvars}}{Logical indicating if standardization was applied.}
#'       \item{\code{bwselect}}{Bandwidth selection strategy used.}
#'       \item{\code{method}}{Bandwidth estimation method.}
#'       \item{\code{vce}}{Variance estimation method.}
#'       \item{\code{scaleregul}}{Scaling factor for regularization.}
#'       \item{\code{scalebiascrct}}{Scaling factor for bias correction.}
#'       \item{\code{N}}{Total sample size \eqn{N}.}
#'     }
#'   }
#' }
#' @seealso \code{\link{rd2d}}, \code{\link{print.rdbw2d}}, \code{\link{summary.rdbw2d}}
#' 
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @references
#' \itemize{
#' \item{\href{https://mdcattaneo.github.io/papers/Cattaneo-Titiunik-Yu_2025_BoundaryRD.pdf}{Cattaneo, M. D., Titiunik, R., Yu, R. R. (2025).}
#' Estimation and Inference in Boundary Discontinuity Designs}
#' }
#' 
#' @examples
#' # Simulated example
#' set.seed(123)
#' n <- 5000
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' t <- as.numeric(X1 > 0)
#' Y <- 3 + 2 * X1 + 1.5 * X2 + t + rnorm(n)
#' X <- cbind(X1, X2)
#' b <- matrix(c(0, 0, 0, 1), ncol = 2)
#'
#' # MSE optimal bandwidth for rd2d
#' bws <- rdbw2d(Y, X, t, b)
#'
#' # View the bandwidth selection results
#' print(bws)
#' summary(bws)
#' @export

rdbw2d <- function(Y, X, t, b, p = 1, deriv = c(0,0), tangvec = NULL,
                   kernel = c("tri","triangular","epa","epanechnikov","uni","uniform","gau","gaussian"),
                   kernel_type = c("prod","rad"),
                   bwselect = c("mserd", "imserd", "msetwo", "imsetwo"),
                   method = c("dpi", "rot"), vce = c("hc1","hc0","hc2","hc3"),
                   bwcheck = 20, masspoints = c("check","adjust","off"),
                   C = NULL, scaleregul = 1, scalebiascrct = 1,
                   stdvars = FALSE){

  # Input error handling

  bwselect <- match.arg(bwselect)
  kernel <- match.arg(kernel)
  kernel_type <- match.arg(kernel_type)
  method <- match.arg(method)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  verbose <- FALSE

  d <- t # renaming the variable
  
  # Check Errors

  exit=0

  if (length(Y) != length(d) || length(Y) != nrow(X)) {
    print("Y, d, and rows of X must have the same length")
    exit <- 1
  }

  if (ncol(X) != 2) {
    print("X must have exactly 2 columns")
    exit <- 1
  }

  if (!(is.logical(d) || all(d %in% c(0, 1)))) {
    print("d must be a logical vector or a numeric vector containing only 0 and 1")
    exit <- 1
  }

  if (kernel!="gau" & kernel!="gaussian" & kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    print("kernel incorrectly specified")
    exit <- 1
  }

  if (!is.numeric(deriv) | length(deriv) != 2) {
    print("deriv must be a numeric vector of length 2")
    exit <- 1
  } else if (sum(deriv) > p) {
    print("Sum of deriv components must be less than or equal to polynomial order p")
    exit <- 1
  }

  if (!is.null(tangvec)) {
    if (!(is.matrix(tangvec) || is.data.frame(tangvec)) ||
        nrow(tangvec) != nrow(b) || ncol(tangvec) != 2) {
      print("tangvec must be a matrix or data frame with the same number of rows as b and exactly 2 columns")
      exit <- 1
    }
  }

  if (!is.null(C) && !(vce %in% c("hc0", "hc1"))) {
    warning("When C is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }

  if (exit>0) stop()

  # Data Cleaning

  dat <- cbind(X[,1], X[,2], Y, d)
  dat <- as.data.frame(dat)
  colnames(dat) <- c("x.1", "x.2", "y", "d")
  eval <- as.data.frame(b)
  colnames(eval) <- c("x.1", "x.2")
  neval <- dim(eval)[1]
  na.ok <- complete.cases(dat$x.1) & complete.cases(dat$x.2) & complete.cases(dat$y) & complete.cases(dat$d)
  dat <- dat[na.ok,]
  N <- dim(dat)[1]
  N.0 <- dim(dat[dat$d == 0,])[1]
  N.1 <- dim(dat[dat$d == 1,])[1]

  if (is.null(p))         p <- 1
  kernel   <- tolower(kernel)

  e_deriv <- matrix(0, nrow = neval, ncol = factorial(p+2)/(factorial(p) * factorial(2)))
  deriv.sum <- deriv[1] + deriv[2]
  if (deriv.sum >= 1){
    e_deriv[,(factorial(deriv.sum+1)/(factorial(deriv.sum-1)* 2)) + deriv[2] + 1] <- 1
  } else {
    e_deriv[,1] <- 1
  }


  if (!is.null(tangvec)){
    warning("Tangvec provided. Ignore option deriv.")
    e_deriv <- matrix(0, nrow = neval, ncol = factorial(p+2)/(factorial(p) * factorial(2)))
    e_deriv[,2] <- tangvec[,1]
    e_deriv[,3] <- tangvec[,2]
    deriv <- c(1,0) # standardization for latter codes
    deriv.sum <- 1
  }

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  # Standardize data if necessary

  if (stdvars){
    sd.1 <- sd(dat$x.1)
    sd.2 <- sd(dat$x.2)
    dat$x.1 <- dat$x.1 / sd.1
    dat$x.2 <- dat$x.2 / sd.2
    eval$x.1 <- eval$x.1 / sd.1
    eval$x.2 <- eval$x.2 / sd.2
  } else {
    sd.1 <- 1
    sd.2 <- 1
  }

  # Store variance and bias constants for IMSE

  bconst <- rep(NA, neval)
  vconst <- rep(NA, neval)

  # Check for mass points

  M <- N; M.0 <- N.0; M.1 <- N.1
  if (masspoints == "check" | masspoints == "adjust"){
    unique.const <- rd2d_unique(dat)
    unique <- unique.const$unique
    M.0 <- dim(unique[unique$d == 0,])[1]
    M.1 <- dim(unique[unique$d == 1,])[1]
    M <- M.0 + M.1
    mass <- 1 - M / N
    if (mass >= 0.2){
      warning("Mass points detected in the running variables.")
      if (masspoints == "check") warning("Try using option masspoints=adjust.")
      if (is.null(bwcheck) & (masspoints == "check" | masspoints == "adjust")) bwcheck <- 50 + p + 1
    }
  }
  # if (!is.null(C)){
  #   unique.0 <- unique(C[d == FALSE])
  #   unique.1 <- unique(C[d == TRUE])
  #   M.0 <- length(unique.0)
  #   M.1 <- length(unique.1)
  #   M <- M.0 + M.1
  # }

  # Rule of thumb bandwidth selection

  dn <- rdbw2d_rot(dat,kernel.type, M)

  # Loop over points of evaluations
  results <- data.frame(matrix(NA, ncol = 16, nrow = neval))
  colnames(results) <- c('b1','b2','h01', 'h02', 'h11', 'h12', 'Nh.0', 'Nh.1',
                         'bias.0', 'bias.1', 'var.0', 'var.1','reg.bias.0','reg.bias.1','reg.var.0','reg.var.1')

  for (i in 1:neval){

    ev <- eval[i,]
    vec <- e_deriv[i,]

    # Center data

    dat.centered <- dat[,c("x.1", "x.2", "y", "d")]
    dat.centered$x.1 <- dat.centered$x.1 - ev$x.1
    dat.centered$x.2 <- dat.centered$x.2 - ev$x.2
    dat.centered$dist <- sqrt(dat.centered$x.1^2 + dat.centered$x.2^2)

    if (masspoints == "adjust"){
      unique.centered <- unique
      unique.centered$x.1 <- unique.centered$x.1 - ev$x.1
      unique.centered$x.2 <- unique.centered$x.2 - ev$x.2
      unique.centered$dist <- sqrt(unique.centered$x.1^2 + unique.centered$x.2^2)
    }

    # Weights

    dn.0 <- dn; dn.1 <- dn

    if (!is.null(bwcheck)) { # Bandwidth restrictions
      if (masspoints == "adjust"){
        sorted.0 <- sort(unique.centered[unique.centered$d == FALSE,]$dist)
        sorted.1 <- sort(unique.centered[unique.centered$d == TRUE,]$dist)
        bw.min.0   <- sorted.0[bwcheck]
        bw.min.1   <- sorted.1[bwcheck]
        bw.max.0   <- sorted.0[length(sorted.0)]
        bw.max.1   <- sorted.1[length(sorted.1)]
      } else{
        sorted.0   <- sort(dat.centered[dat.centered$d == FALSE,]$dist)
        sorted.1   <- sort(dat.centered[dat.centered$d == TRUE,]$dist)
        bw.min.0   <- sorted.0[bwcheck]
        bw.min.1   <- sorted.1[bwcheck]
        bw.max.0   <- sorted.0[length(sorted.0)]
        bw.max.1   <- sorted.1[length(sorted.1)]
      }
      dn.0     <- max(dn, bw.min.0)
      dn.1     <- max(dn, bw.min.1)
      dn.0     <- min(dn, bw.max.0)
      dn.1     <- min(dn, bw.max.1)
    }

    if (kernel_type == "prod"){
      w.0 <- W.fun(dat.centered[dat.centered$d == FALSE,]$x.1/dn.0, kernel) *
             W.fun(dat.centered[dat.centered$d == FALSE,]$x.2/dn.0, kernel) / c(dn.0^2)
      w.1 <- W.fun(dat.centered[dat.centered$d == TRUE,]$x.1/dn.1, kernel) *
             W.fun(dat.centered[dat.centered$d == TRUE,]$x.2/dn.1, kernel) / c(dn.1^2)
    }
    else{
      w.0   <- W.fun(dat.centered[dat.centered$d == FALSE,]$dist/dn.0, kernel)/c(dn.0^2)
      w.1   <- W.fun(dat.centered[dat.centered$d == TRUE,]$dist/dn.1, kernel)/c(dn.1^2)
    }

    eN.0 <- sum(w.0 > 0)
    eN.1 <- sum(w.1 > 0)

    vec.q.0 <- get_coeff(dat.centered[dat.centered$d == FALSE,], vec, p, dn, kernel, kernel_type)
    vec.q.1 <- get_coeff(dat.centered[dat.centered$d == TRUE,], vec, p, dn, kernel, kernel_type)

    if (verbose) {print("Coefficients for a linear combination of (p+1)-th derivatives"); print(vec.q.0); print(vec.q.1)}

    # Bandwidth for fitting the linear combination of (p+1)-th derivatives using (p+1)-th degree model.

    thrshd.0 <- median(dat.centered[dat.centered$d == FALSE,]$dist)
    thrshd.1 <- median(dat.centered[dat.centered$d == TRUE,]$dist)

    bn.0 <- thrshd.0 # If method is "rot", use half of control data to estimate (p+1)th derivative of control.
    bn.1 <- thrshd.1 # If method is "rot", use half of treated data to estimate (p+1)th derivative of treated.

    if (method == "dpi"){

      bn.const.0 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == FALSE,], p + 1, vec.q.0, dn.0, thrshd.0, NULL, vce, kernel, kernel_type, C[as.logical(dat.centered$d == FALSE)])
      bn.const.1 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == TRUE,], p + 1, vec.q.1, dn.1, thrshd.1, NULL, vce, kernel, kernel_type,C[as.logical(dat.centered$d == TRUE)])

      bn.0 <-  ((2 + 2 * (p+1)) * bn.const.0$V  / ( (2 * (p + 1) + 2 - 2 * (p+1)) * (bn.const.0$B^2 + scaleregul * bn.const.0$Reg.1) ) )^(1/(2 * p + 6))
      bn.1 <-  ((2 + 2 * (p+1)) * bn.const.1$V  / ( (2 * (p + 1) + 2 - 2 * (p+1)) * (bn.const.1$B^2 + scaleregul * bn.const.1$Reg.1) ) )^(1/(2 * p + 6))

      if (verbose) print(paste("bn.0 = ", bn.0, ", bn.1 = ", bn.1, sep = ""))
      if (verbose) {print("Constants for bn.0 and bn.1:"); print(paste("B.0 = ", bn.const.0$B, ", V.0 = ", bn.const.0$V, ", Reg.0 = ", bn.const.0$Reg.1));
        print(paste("B.1 = ", bn.const.1$B, ", V.1 = ", bn.const.1$V, ", Reg.1 = ", bn.const.1$Reg.1))}

      if (!is.null(bwcheck)){ # Bandwidth restrictions
        bn.0 <- max(bn.0, bw.min.0)
        bn.1 <- max(bn.1, bw.min.1)
        bn.0 <- min(bn.0, bw.max.0)
        bn.1 <- min(bn.1, bw.max.1)
      }
    }

    # Bandwidth for estimating (derivatives of) treatment effect using p-th degree model.

    if (bwselect == "mserd" | bwselect == "imserd"){

      hn.const.0 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == FALSE,], p, vec, dn.0, bn.0, thrshd.0, vce, kernel, kernel_type, C[as.logical(dat.centered$d == FALSE)])
      hn.const.1 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == TRUE,], p, vec, dn.1, bn.1, thrshd.1, vce, kernel, kernel_type, C[as.logical(dat.centered$d == TRUE)])

      hn <- ( (2 + 2 * deriv.sum) * (hn.const.0$V   + hn.const.1$V) /
                ( (2 * p + 2 - 2 * deriv.sum) * ( (hn.const.0$B + scalebiascrct * hn.const.0$Reg.2 - hn.const.1$B - scalebiascrct * hn.const.1$Reg.2)^2 +
                                     scaleregul * hn.const.0$Reg.1 + scaleregul * hn.const.1$Reg.1) ) )^(1/(2 * p + 4))

      if (!is.null(bwcheck)) { # Bandwidth restrictions
        hn     <- max(hn, bw.min.0, bw.min.1)
        hn     <- min(hn, max(bw.max.0, bw.max.1))
      }

      hn.0 <- hn.1 <- hn
    }

    if (bwselect == "msetwo" | bwselect == "imsetwo"){

      hn.const.0 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == FALSE,], p, vec, dn.0, bn.0, thrshd.0, vce, kernel, kernel_type, C[as.logical(dat.centered$d == FALSE)])
      hn.const.1 <- rdbw2d_bw_v2(dat.centered[dat.centered$d == TRUE,], p, vec, dn.1, bn.1, thrshd.1, vce, kernel, kernel_type, C[as.logical(dat.centered$d == TRUE)])

      hn.0 <- ( (2 + 2 * deriv.sum) * hn.const.0$V /
                ( (2 * p + 2) * ( (hn.const.0$B + scalebiascrct * hn.const.0$Reg.2)^2 + scaleregul * hn.const.0$Reg.1) ) )^(1/(2 * p + 4))
      hn.1 <- ( (2 + 2 * deriv.sum) * hn.const.1$V /
                ( (2 * p + 2) * ( (hn.const.1$B + scalebiascrct * hn.const.1$Reg.2)^2 + scaleregul * hn.const.1$Reg.1) ) )^(1/(2 * p + 4))

      if (!is.null(bwcheck)) { # Bandwidth restrictions
        hn.0     <- max(hn.0, bw.min.0)
        hn.1     <- max(hn.1, bw.min.1)
        hn.0     <- min(hn.0, bw.max.0)
        hn.1     <- min(hn.1, bw.max.1)
      }
    }

    results[i,c(1:2)] <- c(ev$x.1, ev$x.2)
    results[i,c(3:16)] <- c(hn.0, hn.0, hn.1, hn.1, eN.0, eN.1, hn.const.0$B, hn.const.1$B, hn.const.0$V, hn.const.1$V, hn.const.0$Reg.2, hn.const.1$Reg.2,
                            hn.const.0$Reg.1, hn.const.1$Reg.1)
  }

  if (bwselect == "imserd"){
    V.V <- mean(results$var.0) + mean(results$var.1)
    B.B <- mean( (results$bias.0 + scalebiascrct * results$reg.bias.0 - results$bias.1 - scalebiascrct * results$reg.bias.1)^2 + scaleregul * results$reg.var.0 + scaleregul * results$reg.var.1 )
    hIMSE <- ((2 + 2 * deriv.sum) * V.V / ( (2 * p + 2) * B.B ) )^(1/(2 * p + 4))
    results$h01 <- rep(hIMSE, dim(results)[1])
    results$h02 <- rep(hIMSE, dim(results)[1])
    results$h11 <- rep(hIMSE, dim(results)[1])
    results$h12 <- rep(hIMSE, dim(results)[1])
  }

  if (bwselect == "imsetwo"){
    V.V.0 <- mean(results$var.0)
    V.V.1 <- mean(results$var.1)
    B.B.0 <- mean( (results$bias.0 + scalebiascrct * results$reg.bias.0)^2 + scaleregul * results$reg.var.0)
    B.B.1 <- mean( (results$bias.1 + scalebiascrct * results$reg.bias.1)^2 + scaleregul * results$reg.var.1)
    hIMSE.0 <- ((2 + 2 * deriv.sum) * V.V.0 / ( (2 * p + 2) * B.B.0 ) )^(1/(2 * p + 4))
    hIMSE.1 <- ((2 + 2 * deriv.sum) * V.V.1 / ( (2 * p + 2) * B.B.1 ) )^(1/(2 * p + 4))
    results$h01 <- rep(hIMSE.0, dim(results)[1])
    results$h02 <- rep(hIMSE.0, dim(results)[1])
    results$h11 <- rep(hIMSE.1, dim(results)[1])
    results$h12 <- rep(hIMSE.1, dim(results)[1])
  }

  # Standardization (sd.1 = sd.2 = 1 if stdvar == FALSE)

  results$h01 <- results$h01 * sd.1
  results$h02 <- results$h02 * sd.2
  results$h11 <- results$h11 * sd.1
  results$h12 <- results$h12 * sd.2

  # Outputs

  bws <- results[,c("h01", "h02", "h11", "h12")]
  bws <- cbind(eval, bws)
  colnames(bws) <- c("b1","b2","h01", "h02", "h11", "h12")

  clustered <- !is.null(C)

  out        <- list(bws = bws, mseconsts = results,
                     opt = list(N=N, N.0 = N.0, N.1 = N.1,M.0 = M.0,
                                           M.1 = M.1, neval=neval, p=p, deriv=deriv, tangvec = tangvec,
                                           kernel=kernel.type, kernel_type = kernel_type,
                                           bwselect=bwselect, method = method, bwcheck = bwcheck,
                                           stdvars = stdvars, C= C, clustered = clustered,
                                           vce = vce, masspoints = masspoints,
                                           scaleregul = scaleregul, scalebiascrct = scalebiascrct))
  out$call   <- match.call()
  class(out) <- "rdbw2d"
  return(out)
}

################################################################################
#' Print Method for Bandwidth Selection for 2D Local Polynomial RD Design
#'
#' @description The print method for bandwidth selection for 2D local polynomial RD design
#'
#' @param x Class \code{rdbw2d} objects, obtained by calling \code{\link{rdbw2d}}.
#' @param ... Additional arguments passed to the method (currently ignored).
#'
#' @return No return value, called to print \code{\link{rdbw2d}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso \code{\link{rdbw2d}} for bandwidth selection for 2D local polynomial RD design
#'
#' Supported methods: \code{\link{print.rdbw2d}}, \code{\link{summary.rdbw2d}}.
#'
#' @export
#'

print.rdbw2d <- function(x,...){
  cat("Call: rdbw2d\n\n")

  # Format and print the vector as "(x, y)"
  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  cat(sprintf("BW type.               %s\n", paste(x$opt$bwselect, x$opt$method, sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), x$opt$kernel_type, sep = "-")))
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  cat(sprintf("Standardization        %s\n", ifelse(x$opt$stdvars, "on", "off")))
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Estimand (deriv)       %-10d   %-10d\n", x$opt$deriv[1], x$opt$deriv[2]))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  if (x$opt$masspoints == "check" | x$opt$masspoints == "adjust") {
    cat(sprintf("Unique Obs.            %-10d   %-10d\n", x$opt$M.0, x$opt$M.1))
  }
  cat("\n")
}

################################################################################
#' Summary Method for Bandwidth Selection for 2D Local Polynomial RD Design
#'
#' @description
#' Summary method for objects of class \code{rdbw2d}, displaying bandwidth selection results for 2D local polynomial regression discontinuity designs.
#'
#' @param object An object of class \code{rdbw2d}, typically returned by \code{\link{rdbw2d}}.
##' @param ... Optional arguments. Supported options include:
#'   \itemize{
#'     \item \code{subset}: Integer vector of indices of evaluation points to display. Defaults to all evaluation points.
#'     \item \code{sep}: Integer. Controls spacing in the output. Default is \code{8}.
#'   }
#'   
#' @return No return value. Called for its side effects of printing a formatted summary of \code{\link{rdbw2d}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso \code{\link{rdbw2d}} for bandwidth selection in 2D local polynomial RD design.
#'
#' Supported methods: \code{\link{print.rdbw2d}}, \code{\link{summary.rdbw2d}}.
#'
#' @export

summary.rdbw2d <- function(object, ...) {

  x <- object
  
  args <- list(...)
  
  if (is.null(args[['subset']])) {
    subset <- NULL
  } else {
    subset <- args[['subset']]
  }
  
  if (is.null(args[['sep']])) {
    sep <- 8
  } else {
    sep <- args[['sep']]
  }

  cat("Call: rdbw2d\n\n")

  # Format and print the vector as "(x, y)"
  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  cat(sprintf("BW type.               %s\n", paste(x$opt$bwselect, x$opt$method, sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), x$opt$kernel_type, sep = "-")))
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  cat(sprintf("Standardization        %s\n", ifelse(x$opt$stdvars, "on", "off")))
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Estimand (deriv)       %-10d   %-10d\n", x$opt$deriv[1], x$opt$deriv[2]))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  if (x$opt$masspoints == "check" | x$opt$masspoints == "adjust") {
    cat(sprintf("Unique Obs.            %-10d   %-10d\n", x$opt$M.0, x$opt$M.1))
  }
  cat("\n")

  # Define column headers and widths

  cat("Bandwidth Selection","\n")

  headers <- c("ID", "b1", "b2", "h01", "h02", "h11", "h12")
  col_widths <- c(4, sep, sep, sep, sep, sep, sep)

  cat(strrep("=", sum(col_widths)), "\n")
  group_headers <- c(
    formatC("        Bdy Points", width = col_widths[1] + col_widths[2] + col_widths[3], format = "s", flag = "-"),
    formatC("      BW Control", width = col_widths[4] + col_widths[5], format = "s", flag = "-"),
    formatC("    BW Treatment", width = col_widths[6] + col_widths[7], format = "s", flag = "-")
  )
  cat(paste(group_headers, collapse = ""), "\n")

  # Format and print header row
  formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
  cat(paste(formatted_headers, collapse = ""), "\n")
  cat(strrep("=", sum(col_widths)), "\n")

  neval <- nrow(x$bws)
  if (is.null(subset)){
    subset <- seq_len(neval)
  } else{
    # input error handling
    if (!all(subset %in% seq_len(neval))) {
      warning("Invalid subset provided. Resetting to default: 1:neval")
      subset <- seq_len(neval)
    }
  }

  # Print each row of bandwidth estimates
  for (j in 1:nrow(x$bws)) {
    index <- formatC(j, width = col_widths[1], format = "d")
    bdy1 <- ifelse(is.na(x$bws[j, "b1"]),
                   formatC("NA", width = col_widths[2], format = "s"),
                   formatC(x$bws[j, "b1"], format = "f", digits = 3, width = col_widths[2]))
    bdy2 <- ifelse(is.na(x$bws[j, "b2"]),
                   formatC("NA", width = col_widths[3], format = "s"),
                   formatC(x$bws[j, "b2"], format = "f", digits = 3, width = col_widths[3]))
    control1 <- ifelse(is.na(x$bws[j, "h01"]),
                       formatC("NA", width = col_widths[4], format = "s"),
                       formatC(x$bws[j, "h01"], format = "f", digits = 3, width = col_widths[4]))
    control2 <- ifelse(is.na(x$bws[j, "h02"]),
                       formatC("NA", width = col_widths[5], format = "s"),
                       formatC(x$bws[j, "h02"], format = "f", digits = 3, width = col_widths[5]))
    treatment1 <- ifelse(is.na(x$bws[j, "h11"]),
                         formatC("NA", width = col_widths[6], format = "s"),
                         formatC(x$bws[j, "h11"], format = "f", digits = 3, width = col_widths[6]))
    treatment2 <- ifelse(is.na(x$bws[j, "h12"]),
                         formatC("NA", width = col_widths[7], format = "s"),
                         formatC(x$bws[j, "h12"], format = "f", digits = 3, width = col_widths[7]))

    # Combine formatted values and print the row
    row_vals <- c(index, bdy1, bdy2, control1, control2, treatment1, treatment2)
    if (j %in% subset) cat(paste(row_vals, collapse = ""), "\n")
  }

  # Print closing separator line
  cat(strrep("=", sum(col_widths)), "\n")
}
