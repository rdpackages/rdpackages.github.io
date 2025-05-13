#' MSE bandwidth selection for geometrical RD design
#'
#' @title Local Polynomial RD Estimation on Distance-Based Running Variables
#'
#' @description
#' \code{rd2d.dist} implements distance-based local polynomial boundary regression discontinuity (RD) point estimators with robust bias-corrected pointwise confidence intervals and 
#' uniform confidence bands, developed in Cattaneo, Titiunik and Yu (2025). For robust bias-correction, see Calonico, Cattaneo and Titiunik (2014).
#' 
#' Companion commands are: \code{rdbw2d.dist} for data-driven bandwidth selection.
#' 
#' For other packages of RD designs, visit
#' <https://rdpackages.github.io/>
#' @param Y Dependent variable; a numeric vector of length \eqn{N}, where \eqn{N} is the sample size.
#' @param D Distance-based scores \eqn{\mathbf{D}_i=(\mathbf{D}_{i}(\mathbf{b}_1),\cdots,\mathbf{D}_{i}(\mathbf{b}_J))}; dimension is \eqn{N \times J} where \eqn{N} = sample size and \eqn{J} = number of cutoffs; non-negative values means data point in treatment group and negative values means data point in control group.
#' @param h Bandwidth(s); if \eqn{c=h} then same bandwidth is used for both groups; if a matrix of size \eqn{J \times 2} is provided, each row contains \eqn{(h_{\text{control}}, h_{\text{tr}})} for the evaluation point; if not specified, bandwidths are selected via \code{rdbw2d.dist()}.
#' @param b Optional evaluation points; a matrix or data frame specifying boundary points \eqn{\mathbf{b}_j = (b_{1j}, b_{2j})}, dimension \eqn{J \times 2}.
#' @param p Polynomial order for point estimation. Default is \code{p = 1}.
#' @param q Polynomial order for bias-corrected estimation. Must satisfy \eqn{q \geq p}. Default is \code{q = p + 1}.
#' @param kink Logical; whether to apply kink adjustment. Options: \code{"on"} (default) or \code{"off"}.
#' @param kernel Kernel function to use. Options are \code{"unif"}, \code{"uniform"} (uniform), \code{"triag"}, \code{"triangular"} (triangular, default), and \code{"epan"}, \code{"epanechnikov"} (Epanechnikov).
#' @param level Nominal confidence level for intervals/bands, between 0 and 100 (default is 95).
#' @param cbands Logical. If \code{TRUE}, also compute uniform confidence bands (default is \code{FALSE}).
#' @param side Type of confidence interval. Options: \code{"two"} (two-sided, default), \code{"left"} (left tail), or \code{"right"} (right tail).
#' @param repp Number of bootstrap repetitions used for critical value simulation. Default is \code{1000}.
#' @param bwselect Bandwidth selection strategy. Options: 
#' \itemize{
#' \item \code{"mserd"}. One common MSE-optimal bandwidth selector for the boundary RD treatment effect estimator for each evaluation point (default). 
#' \item \code{"imserd"}. IMSE-optimal bandwidth selector for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"msetwo"}. Two different MSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator for each evaluation point. 
#' \item \code{"imsetwo"}. Two IMSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"user provided"}. User-provided bandwidths. If \code{h} is not \code{NULL}, then \code{bwselect} is overwritten to \code{"user provided"}.
#' }
#' @param vce Variance-covariance estimator for standard errors.
#' Options:
#' \describe{
#'   \item{\code{"hc0"}}{Heteroskedasticity-robust variance estimator without small sample adjustment (White robust).}
#'   \item{\code{"hc1"}}{Heteroskedasticity-robust variance estimator with degrees-of-freedom correction. (default)}
#'   \item{\code{"hc2"}}{Heteroskedasticity-robust variance estimator using leverage adjustments.}
#'   \item{\code{"hc3"}}{More conservative heteroskedasticity-robust variance estimator (similar to jackknife correction).}
#' }
#' @param rbc Logical. Whether to apply robust bias correction. Options: \code{"on"} (default) or \code{"off"}. When \code{kink = off}, turn on \code{rbc} means setting \code{q} to \code{p + 1}. 
#' When \code{kink = on}, turn on \code{rbc} means shrinking the bandwidth selector to be proportional to \eqn{N^{-1/3}}.
#' @param bwcheck If a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least \code{bwcheck} unique observations are used. Default is \code{50 + p + 1}.
#' @param masspoints Strategy for handling mass points in the running variable.
#' Options:
#' \describe{
#'   \item{\code{"check"}}{(default) Check for repeated values and adjust inference if needed.}
#'   \item{\code{"adjust"}}{Adjust bandwidths to guarantee a sufficient number of unique support points.}
#'   \item{\code{"off"}}{Ignore mass points completely.}
#' }
#' @param C Cluster ID variable used for cluster-robust variance estimation with degrees-of-freedom weights.Default is \code{C = NULL}.
#' @param scaleregul Scaling factor for the regularization term in bandwidth selection. Default is \code{1}.
#' @param cqt Constant controlling subsample fraction for initial bias estimation. Default is \code{0.5}.
#'
#' @return An object of class \code{"rd2d.dist"}, a list containing:
#' \describe{
#'   \item{\code{main}}{Data frame of point estimates, standard errors, confidence intervals, and bandwidths:
#'     \describe{
#'       \item{\code{b1}}{First coordinate of the evaluation point.}
#'       \item{\code{b2}}{Second coordinate of the evaluation point.}
#'       \item{\code{Est.p}}{Point estimate \eqn{\widehat{\tau}_{\text{dist},p}(\mathbf{b})} with polynomial order \eqn{p}.}
#'       \item{\code{Var.p}}{Variance of \eqn{\widehat{\tau}_{\text{dist},p}(\mathbf{b})}.}
#'       \item{\code{Est.q}}{Bias-corrected estimate \eqn{\widehat{\tau}_{\text{dist},q}(\mathbf{b})} with polynomial order \eqn{q}.}
#'       \item{\code{Var.q}}{Variance of \eqn{\widehat{\tau}_{\text{dist},q}(\mathbf{b})}.}
#'       \item{\code{pvalue}}{Two-sided p-value based on \eqn{T_{\text{dist},q}(\mathbf{b})}.}
#'       \item{\code{CI.lower}}{Lower bound of confidence interval.}
#'       \item{\code{CI.upper}}{Upper bound of confidence interval.}
#'       \item{\code{CB.lower}}{Lower bound of uniform confidence band (if \code{cbands=TRUE}).}
#'       \item{\code{CB.upper}}{Upper bound of uniform confidence band (if \code{cbands=TRUE}).}
#'       \item{\code{h0}}{Bandwidth used for control group (\eqn{D_i(\mathbf{b}) < 0}).}
#'       \item{\code{h1}}{Bandwidth used for treatment group (\eqn{D_i(\mathbf{b}) \geq 0}).}
#'       \item{\code{Nh0}}{Effective sample size for control group.}
#'       \item{\code{Nh1}}{Effective sample size for treatment group.}
#'     }
#'   }
#'   \item{\code{main.A0}}{Summary table for the control group only.}
#'   \item{\code{main.A1}}{Summary table for the treatment group only.}
#'   \item{\code{tau.hat}}{Vector of point estimates \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'   \item{\code{se.hat}}{Standard errors corresponding to \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'   \item{\code{cb}}{Confidence intervals and uniform bands.}
#'   \item{\code{cov.us}}{Covariance matrix used to construct uniform bands.}
#'   \item{\code{opt}}{A list of estimation options (e.g., \code{p}, \code{q}, \code{kernel}, \code{level}, etc.) and internal variables such as sample size \eqn{N}.}
#' }
#'
#' @seealso \code{\link{rdbw2d.dist}}, \code{\link{rd2d}}, \code{\link{print.rd2d.dist}}, \code{\link{summary.rd2d.dist}}
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
#' set.seed(123)
#' n <- 5000
#'
#' # Generate running variables x1 and x2
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#'
#' # Define treatment assignment: treated if x1 >= 0
#' d <- as.numeric(x1 >= 0)
#'
#' # Generate outcome variable y with some treatment effect
#' y <- 3 + 2 * x1 + 1.5 * x2 + 1.5 * d + rnorm(n, sd = 0.5)
#'
#' # Define evaluation points (e.g., at the origin and another point)
#' eval <- data.frame(x.1 = c(0, 0), x.2 = c(0, 1))
#'
#' # Compute Euclidean distances to evaluation points
#' dist.a <- sqrt((x1 - eval$x.1[1])^2 + (x2 - eval$x.2[1])^2)
#' dist.b <- sqrt((x1 - eval$x.1[2])^2 + (x2 - eval$x.2[2])^2)
#'
#' # Combine distances into a matrix
#' D <- as.data.frame(cbind(dist.a, dist.b))
#'
#' # Assign positive distances for treatment group, negative for control
#' d_expanded <- matrix(rep(2 * d - 1, times = ncol(D)), nrow = nrow(D), ncol = ncol(D))
#' D <- D * d_expanded
#'
#' # Run the rd2d.dist function
#' result <- rd2d.dist(y, D, b = eval)
#'
#' # View the estimation results
#' print(result)
#' summary(result)
#' @export

rd2d.dist <- function(Y, D, h = NULL, b = NULL, p = 1, q = 2, kink = c("off", "on"),
                      kernel = c("tri","triangular", "epa","epanechnikov","uni","uniform","gau","gaussian"),
                      level = 95, cbands = TRUE, side = c("two", "left", "right"), repp = 1000,
                      bwselect = c("mserd", "imserd", "msetwo", "imsetwo", "user provided"),
                      vce = c("hc1","hc0","hc2","hc3"), rbc = c("on", "off"),
                      bwcheck = 50 + p + 1, masspoints = c("check","adjust","off"),
                      C = NULL, scaleregul = 1, cqt = 0.5){

  # Input error handling

  bwselect <- match.arg(bwselect)
  kernel <- match.arg(kernel)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  side <- match.arg(side)
  kink <- match.arg(kink)
  rbc <- match.arg(rbc)

  # Check Errors

  exit=0

  if (length(Y) != nrow(D)) {
    print("Y and rows of D must have the same length")
    exit <- 1
  }

  if (kernel!="gau" & kernel!="gaussian" & kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    print("kernel incorrectly specified")
    exit <- 1
  }

  if (!is.null(b)){
    if (nrow(b) != ncol(D) || ncol(b) != 2){
      print("b must have 2 columns and the same number of rows as D's number of columns")
      exit <- 1
    }
  }

  # level must be numeric in (0, 100)
  if (!is.numeric(level) || level <= 0 || level >= 100) {
    print("level must be a numeric value between 0 and 100")
    exit <- 1
  }

  # repp must be a positive integer
  if (!is.numeric(repp) || repp < 1 || repp != as.integer(repp)) {
    print("repp must be a positive integer")
    exit <- 1
  }

  # h must be either a positive scalar or a matrix/data.frame with same rows as b and 4 columns
  if (!is.null(h)) {
    if (length(h) == 1) {
      if (!is.numeric(h) || h <= 0) {
        print("If h is a scalar, it must be a positive numeric value")
        exit <- 1
      }
    } else if (!(is.matrix(h) || is.data.frame(h)) ||
               nrow(h) != ncol(D) || ncol(h) != 2) {
      print("If h is not a scalar, it must be a matrix or data frame with the same number of rows as b and 2 columns")
      exit <- 1
    }
  }

  if (is.null(h) & bwselect == "user provided"){
    exit <- 1
    print("Please provide bandwidths.")
  }

  if (!is.null(C) && !(vce %in% c("hc0", "hc1"))) {
    warning("When C is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }

  if (exit>0) stop()

  ############################ Data Preparation ################################

  neval <- ncol(D)

  if (!is.null(b)){
    eval <- as.data.frame(b)
  } else {
    eval <- matrix(NA, nrow = neval, ncol = 2)
    eval <- as.data.frame(eval)
  }

  colnames(eval) <- c("x.1", "x.2")

  # D <- as.data.frame(D)
  na.ok <- complete.cases(D)
  D <- D[na.ok,,drop = FALSE]
  Y <- Y[na.ok]
  C <- C[na.ok]
  d <- (D[,1] >= 0)

  N <- length(Y)
  N.1 <- sum(d)
  N.0 <- N - N.1

  if (is.null(p))         p <- 1
  kernel   <- tolower(kernel)

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  ################################ Bandwidth ###################################

  if (is.null(h)){
    bws <- rdbw2d.dist(Y = Y, D = D, b = b, p = p, kink = kink, kernel = kernel, bwselect = bwselect,
                  vce = vce, bwcheck = bwcheck, masspoints = masspoints,
                  C = C, scaleregul = scaleregul, cqt = cqt)
    bws <- bws$bws
    hgrid <- bws[,3]
    hgrid.1 <- bws[,4]
  } else {
    bwselect <- "user provided"
    # standardize bandwidth
    if (length(h) == 1){
      hgrid <- rep(h, neval)
      hgrid.1 <- rep(h, neval)
    } else {
      hgrid <- h[,1]
      hgrid.1 <- h[,2]
    }
  }

  hfull <- cbind(hgrid, hgrid.1)

  ###################### Point estimation and inference ########################

  # kink adjustment

  # if (kink =="on"){
  #   if (bwselect == "mserd" | bwselect == "imserd"){
  #     hfull <- hfull* N^(-1/4) / N^(-1/(2 * p + 4))
  #   }
  #   if (bwselect == "msetwo" | bwselect == "imsetwo"){
  #     hfull[,1] <- hfull[,1] * N.0^(-1/4) / N.0^(-1/(2 * p + 4))
  #     hfull[,2] <- hfull[,2] * N.1^(-1/4) / N.1^(-1/(2 * p + 4))
  #   }
  #   # if user provides bandwidth, then do not adjust
  # }

  distfit.p <- rd2d_dist_fit(Y = Y, D = D, h = hfull, p = p, b = b, kernel = kernel,
                              vce = vce, bwcheck = bwcheck, masspoints = masspoints, C = C, cbands = FALSE)
  estimate.p <- distfit.p$Estimate
  tau.hat.p <- estimate.p$mu1 - estimate.p$mu0
  se.hat.p <- sqrt(estimate.p$se0^2 + estimate.p$se1^2)
  h0.p <- estimate.p$h0
  h1.p <- estimate.p$h1
  eN0.p <- estimate.p$N0
  eN1.p <- estimate.p$N1

  M.vec <- distfit.p$M.vec
  M.0.vec <- distfit.p$M.0.vec
  M.1.vec <- distfit.p$M.1.vec

  # robust bias correction

  hfull.rbc <- hfull
  if (kink == "on"){
    if (q > p){
      q <- p
      print("q is taken to be p with kink on.")
    }
    if (rbc == "on"){
      if (bwselect == "mserd" | bwselect == "imserd"){
        hfull.rbc <- hfull* M.vec^(-1/3) / M.vec^(-1/4)
      }
      if (bwselect == "msetwo" | bwselect == "imsetwo"){
        hfull.rbc[,1] <- hfull[,1] * M.0.vec^(-1/3) / M.0.vec^(-1/4)
        hfull.rbc[,2] <- hfull[,2] * M.1.vec^(-1/3) / M.1.vec^(-1/4)
      }
    }
  } else{
    if (q == p){
      if (rbc == "on"){
        q <- p + 1
        print("Set q = p + 1 for robust bias correction.")
      }
    }
  }

  distfit.q <- rd2d_dist_fit(Y = Y, D = D, h = hfull.rbc, p = q, b = b, kernel = kernel,
                             vce = vce, bwcheck = bwcheck, masspoints = masspoints, C = C, cbands = cbands)
  estimate.q <- distfit.q$Estimate
  tau.hat.q <- estimate.q$mu1 - estimate.q$mu0
  se.hat.q <- sqrt(estimate.q$se0^2 + estimate.q$se1^2)
  h0.q <- estimate.q$h0
  h1.q <- estimate.q$h1
  eN0.q <- estimate.q$N0
  eN1.q <- estimate.q$N1

  zvalues <- tau.hat.q/se.hat.q
  pvalues <- 2 * pnorm(abs(zvalues),lower.tail = FALSE)

  if (side == "two"){
    zval <- qnorm((level + 100)/ 200)
    CI.lower <- tau.hat.q - zval * se.hat.q
    CI.upper <- tau.hat.q + zval * se.hat.q
  }
  if (side == "left"){
    zval <- qnorm(level / 100)
    CI.upper <- tau.hat.q + zval * se.hat.q
    CI.lower <- rep(-Inf, length(CI.upper))
  }
  if (side == "right"){
    zval <- qnorm(level / 100)
    CI.lower <- tau.hat.q - zval * se.hat.q
    CI.upper <- rep(Inf, length(CI.lower))
  }

  # to store matrices

  Indicators.0 <- distfit.q$Indicators.0
  Indicators.1 <- distfit.q$Indicators.1
  inv.designs.0 <- distfit.q$inv.designs.0
  inv.designs.1 <- distfit.q$inv.designs.1
  sig.halfs.0 <- distfit.q$sig.halfs.0
  sig.halfs.1 <- distfit.q$sig.halfs.1
  resd.0 <- distfit.q$resd.0
  resd.1 <- distfit.q$resd.1

  ############################## covariance ####################################

  cov.us = matrix(NA,neval,neval)

  clusters <- unique(C)
  g <- length(clusters)
  k <- q + 1

  if (cbands){
    for (i in 1:neval) {
      for (j in i:neval) {
        cov.0 = NA
        cov.1 = NA

        ind.i.0 <- Indicators.0[[i]]
        ind.j.0 <- Indicators.0[[j]]
        ind.0 <- ind.i.0 * ind.j.0
        ind.i.1 <- Indicators.1[[i]]
        ind.j.1 <- Indicators.1[[j]]
        ind.1 <- ind.i.1 * ind.j.1
        slice.i.0 <- as.logical(ind.0[ind.i.0])
        slice.j.0 <- as.logical(ind.0[ind.j.0])
        slice.i.1 <- as.logical(ind.1[ind.i.1])
        slice.j.1 <- as.logical(ind.1[ind.j.1])

        sigma0.half.i <- sig.halfs.0[[i]]
        sigma0.half.j <- sig.halfs.0[[j]]
        sigma1.half.i <- sig.halfs.1[[i]]
        sigma1.half.j <- sig.halfs.1[[j]]

        resd.i.0 <- unname(resd.0[[i]][slice.i.0])
        resd.j.0 <- unname(resd.0[[j]][slice.j.0])
        resd.i.1 <- unname(resd.1[[i]][slice.i.1])
        resd.j.1 <- unname(resd.1[[j]][slice.j.1])

        vec <- rep(0,q+1)
        vec[1] <- 1

        if (is.null(C)){
          cov.0 <- matrix(vec, nrow = 1) %*% t(sigma0.half.i[slice.i.0,,drop = FALSE]) %*% sigma0.half.j[slice.j.0,,drop = FALSE] %*% matrix(vec, ncol=1)
          cov.0 <- cov.0[1,1]
          cov.1 <- matrix(vec, nrow = 1) %*% t(sigma1.half.i[slice.i.1,,drop = FALSE]) %*% sigma1.half.j[slice.j.1,,drop = FALSE] %*% matrix(vec, ncol=1)
          cov.1 <- cov.1[1,1]
        } else {
          M.0 <- matrix(0, k, k)
          M.1 <- matrix(0, k, k)
          for (l in 1:g){
            M.0 <- M.0 + matrix(sigma0.half.i[l,], ncol = 1) %*% matrix(sigma0.half.j[l,], nrow = 1)
            M.1 <- M.1 + matrix(sigma1.half.i[l,], ncol = 1) %*% matrix(sigma1.half.j[l,], nrow = 1)
          }
          cov.0 <- matrix(vec, nrow = 1) %*% M.0 %*% matrix(vec, ncol=1)
          cov.1 <- matrix(vec, nrow = 1) %*% M.1 %*% matrix(vec, ncol=1)
        }
      cov.us[i,j] <- cov.0 + cov.1
      cov.us[j,i] <- cov.0 + cov.1
      }
    }
  }

  cov.hat.q <- NA
  cb.hat.q <- list(CI.l = CI.lower, CI.r = CI.upper, CB.l = rep(NA, length(CI.lower)), CB.r = rep(NA, length(CI.lower)))
  CB.lower <- NA
  CB.upper <- NA
  if (cbands){
    cov.hat.q <- cov.us
    cb.hat.q <- rd2d_cb(tau.hat.q, cov.hat.q, repp, side, level)
    CB.lower <- cb.hat.q$CB.l
    CB.upper <- cb.hat.q$CB.r
  }

  clustered <- !is.null(C)

  ############################### outputs ######################################

  main <- cbind(eval[,1], eval[,2], tau.hat.p, se.hat.p, tau.hat.q, se.hat.q, zvalues, pvalues,
                CI.lower, CI.upper, CB.lower, CB.upper, hfull[,1], hfull[,2],
                hfull.rbc[,1], hfull.rbc[,2], eN0.p, eN1.p)
  main <- as.data.frame(main)
  colnames(main) <- c("b1","b2","Est.p","Var.p","Est.q","Var.q", "z", "P>|z|",
                      "CI.lower","CI.upper","CB.lower", "CB.upper", "h0", "h1",
                      "h0.rbc", "h1.rbc", "Nh0", "Nh1")

  main.A0 <- cbind(eval[,1], eval[,2], estimate.p$mu0, estimate.p$se0, estimate.q$mu0, estimate.q$se0, hfull[,1], hfull.rbc[,1], eN0.p)
  main.A0 <- as.data.frame(main.A0)
  colnames(main.A0) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q","h0", "h0.rbc","Nh0")

  main.A1 <- cbind(eval[,1], eval[,2], estimate.p$mu1, estimate.p$se1, estimate.q$mu1, estimate.q$se1, hfull[,2], hfull.rbc[,2], eN1.p)
  main.A1 <- as.data.frame(main.A1)
  colnames(main.A1) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q","h1", "h1.rbc","Nh1")

  rdmodel <- "rd2d.dist"

  out <- list(main = main, main.A0 = main.A0, main.A1 = main.A1,
              opt=list(b = eval, p = p, q = q, kernel=kernel.type, kink = kink, N=N, N.0 = N.0, rbc = rbc,
                       N.1 = N.1, M = M.vec, M.0 = M.0.vec, M.1 = M.1.vec, neval=neval, bwselect = bwselect,
                       vce = vce, bwcheck = bwcheck, masspoints = masspoints, C = C, clustered = clustered,
                       scaleregul = scaleregul, cqt = cqt, 
                       level = level, repp = repp, side = side,cbands = cbands,
                       h0 = hfull[,1], h1 = hfull[,2], h0.rbc = hfull.rbc[,1], h1.rbc = hfull.rbc[,2],
                       Nh0 = eN0.p, Nh1 = eN1.p),
              tau.hat = tau.hat.p, tau.hat.q = tau.hat.q, se.hat = se.hat.p, cov.us=cov.us, cb = cb.hat.q, pvalues = pvalues, zvalues = zvalues, rdmodel = rdmodel)
  out$call   <- match.call()
  class(out) <- "rd2d.dist"

  return(out)
}

################################################################################
#' Print Method for 2D Local Polynomial RD Estimation (Distance-Based)
#'
#' @description
#' Prints the results of a 2D local polynomial regression discontinuity (RD) estimation using distance-based evaluation, as obtained from \code{\link{rd2d.dist}}.
#'
#' @param x An object of class \code{rd2d.dist}, returned by \code{\link{rd2d.dist}}.
#' @param ... Additional arguments passed to the method (currently ignored).
#'
#' @return
#' No return value. This function is called for its side effects: it prints the \code{\link{rd2d.dist}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso
#' \code{\link{rd2d.dist}} for estimation using distance-based methods in 2D local polynomial RD designs.
#'
#' Supported methods: \code{\link{print.rd2d.dist}}, \code{\link{summary.rd2d.dist}}.
#'
#' @export
#' 
print.rd2d.dist <- function(x,...) {

  cat(paste(x$rdmodel, "\n", sep = ""))
  cat(paste("\n", sep = ""))

  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  cat(sprintf("BW type                %s\n", paste(x$opt$bwselect, "rot", sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), "rad", sep = "-")))
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Estimand (deriv)       %-10d   %-10d\n", 0, 0))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
  cat("\n")
}

################################################################################
#' Summary Method for 2D Local Polynomial RD Estimation (Distance-Based)
#'
#' @description
#' Summarizes estimation and bandwidth results from a 2D local polynomial regression discontinuity (RD) design using distance-based methods, as returned by \code{\link{rd2d.dist}}.
#'
#' @param object An object of class \code{rd2d.dist}, returned by \code{\link{rd2d.dist}}.
#' @param ... Optional arguments. Supported options include:
#'   \itemize{
#'     \item \code{CBuniform}: Logical. If \code{TRUE}, displays uniform confidence bands;
#'       if \code{FALSE} (default), displays pointwise confidence intervals.
#'     \item \code{subset}: Integer vector of indices of evaluation points to display.
#'       Defaults to all evaluation points.
#'     \item \code{output}: Character. Use \code{"main"} to display estimation results,
#'       or \code{"bw"} to display bandwidth information. Default is \code{"main"}.
#'     \item \code{sep}: Integer vector of length three. Controls spacing in the output.
#'       \code{sep[1]} controls spacing for the columns of boundary points, estimation,
#'       z-value, and p-value in the \code{"main"} table.
#'       \code{sep[2]} controls spacing for confidence intervals (or bands) in the \code{"main"} table.
#'       \code{sep[3]} controls spacing for the columns in the \code{"bw"} table.
#'       Default is \code{c(7, 17, 8)}.
#'   }
#'
#' @return No return value. This function is called for its side effects: it prints a formatted summary of \code{\link{rd2d.dist}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso \code{\link{rd2d.dist}} for estimation using distance-based 2D local polynomial RD design.
#'
#' Supported methods: \code{\link{print.rd2d.dist}}, \code{\link{summary.rd2d.dist}}.
#'
#' @export

summary.rd2d.dist <- function(object, ...) {

    x <- object

    args <- list(...)
    
    if (is.null(args[['CBuniform']])) {
      CBuniform <- FALSE
    } else {
      CBuniform <- TRUE
    }
    
    if (is.null(args[['subset']])) {
      subset <- NULL
    } else {
      subset <- args[['subset']]
    }
    
    if (is.null(args[['output']])) {
      output <- "main" 
    } else {
      output <- args[['output']]
    }
    
    if (is.null(args[['sep']])) {
      sep <- c(7,17,8)
    } else {
      sep <- args[['sep']]
    }
    
    cat(paste(x$rdmodel, "\n", sep = ""))
    cat(paste("\n", sep = ""))

    cat(sprintf("Number of Obs.         %d\n", x$opt$N))
    cat(sprintf("BW type                %s\n", paste(x$opt$bwselect, "rot", sep = "-")))
    cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), "rad", sep = "-")))
    cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
    cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
    cat("\n")
    cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
    cat(sprintf("Estimand (deriv)       %-10d   %-10d\n", 0, 0))
    cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
    cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
    cat("\n")

  results <- data.frame(
    b1 = x$opt$b[, 1],
    b2 = x$opt$b[, 2],
    Coef = x$tau.hat,
    Zvalues = x$zvalues,
    Pvalues = x$pvalues,
    CILower = if (CBuniform) x$cb$CB.l else x$cb$CI.l,
    CIUpper = if (CBuniform) x$cb$CB.r else x$cb$CI.r,
    h0 = x$opt$h0,
    h1 = x$opt$h1,
    Nh0 = x$opt$Nh0,
    Nh1 = x$opt$Nh1
  )

  eval.specified <- !all(is.na(x$opt$b)) # TRUE if at least one b is not NA

  neval <- nrow(results)
  if (is.null(subset)) {
    subset <- seq_len(neval)
  } else {
    if (!all(subset %in% seq_len(neval))) {
      warning("Invalid subset provided. Resetting to default: 1:neval")
      subset <- seq_len(neval)
    }
  }

  if (output == "main") {
    if (eval.specified) {
      headers <- c("ID", "b1", "b2", "Est.", "z", "P > |z|", sprintf("%d%% CI", x$opt$level))
      col_widths <- c(4, sep[1], sep[1], sep[1], sep[1], sep[1], sep[2])
    } else {
      headers <- c("ID", "Est.", "z", "P > |z|", sprintf("%d%% CI", x$opt$level))
      col_widths <- c(4, sep[1], sep[1], sep[1], sep[2])
    }
    if (CBuniform) {
      headers[length(headers)] <- sprintf("%d%% Unif. CB", x$opt$level)
    }

    cat(strrep("=", sum(col_widths) + 2 * (length(headers) - 1)), "\n")
    formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
    cat(paste(formatted_headers, collapse = "  "), "\n")
    cat(strrep("=", sum(col_widths) + 2 * (length(headers) - 1)), "\n")

    for (i in seq_len(neval)) {
      if (i %in% subset) {
        if (eval.specified) {
          row_vals <- c(
            formatC(i, width = col_widths[1], format = "d"),
            formatC(results$b1[i], format = "f", digits = 3, width = col_widths[2]),
            formatC(results$b2[i], format = "f", digits = 3, width = col_widths[3]),
            formatC(results$Coef[i], format = "f", digits = 4, width = col_widths[4]),
            formatC(results$Zvalues[i], format = "f", digits = 4, width = col_widths[5]),
            formatC(results$Pvalues[i], format = "f", digits = 4, width = col_widths[6]),
            formatC(
              paste0("[", formatC(results$CILower[i], format = "f", digits = 4),
                     ", ", formatC(results$CIUpper[i], format = "f", digits = 4), "]"),
              width = col_widths[7], format = "s"
            )
          )
        } else {
          row_vals <- c(
            formatC(i, width = col_widths[1], format = "d"),
            formatC(results$Coef[i], format = "f", digits = 4, width = col_widths[2]),
            formatC(results$Zvalues[i], format = "f", digits = 4, width = col_widths[3]),
            formatC(results$Pvalues[i], format = "f", digits = 4, width = col_widths[4]),
            formatC(
              paste0("[", formatC(results$CILower[i], format = "f", digits = 4),
                     ", ", formatC(results$CIUpper[i], format = "f", digits = 4), "]"),
              width = col_widths[5], format = "s"
            )
          )
        }
        cat(paste(row_vals, collapse = "  "), "\n")
      }
    }
    cat(strrep("=", sum(col_widths) + 2 * (length(headers) - 1)), "\n")

  } else if (output == "bw") {
    if (eval.specified) {
      headers <- c("ID", "b1", "b2", "h0", "h1", "Nh0", "Nh1")
      col_widths <- c(4, sep[3], sep[3], sep[3], sep[3], sep[3], sep[3], sep[3], sep[3])
    } else {
      headers <- c("ID", "h0", "h1", "Nh0", "Nh1")
      col_widths <- c(4, sep[3], sep[3], sep[3], sep[3], sep[3], sep[3])
    }

    cat(strrep("=", sum(col_widths)), "\n")
    if (eval.specified) {
      group_headers <- c(
        formatC("        Bdy Points", width = col_widths[1] + col_widths[2] + col_widths[3], format = "s", flag = "-"),
        formatC("      Bandwidths", width = col_widths[4] + col_widths[5], format = "s", flag = "-"),
        formatC("      Eff. N", width = col_widths[6] + col_widths[7], format = "s", flag = "-")
      )
    } else {
      group_headers <- c(
        formatC("  ", width = col_widths[1], format = "s", flag = "-"),
        formatC("   Bandwidths", width = col_widths[2] + col_widths[3], format = "s", flag = "-"),
        formatC("      Eff. N", width = col_widths[4] + col_widths[5], format = "s", flag = "-")
      )
    }
    cat(paste(group_headers, collapse = ""), "\n")

    formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
    cat(paste(formatted_headers, collapse = ""), "\n")
    cat(strrep("=", sum(col_widths)), "\n")

    for (j in seq_len(neval)) {
      if (j %in% subset) {
        if (eval.specified) {
          row_vals <- c(
            formatC(j, width = col_widths[1], format = "d"),
            formatC(results$b1[j], format = "f", digits = 3, width = col_widths[2]),
            formatC(results$b2[j], format = "f", digits = 3, width = col_widths[3]),
            formatC(results$h0[j], format = "f", digits = 3, width = col_widths[4]),
            formatC(results$h1[j], format = "f", digits = 3, width = col_widths[5]),
            formatC(results$Nh0[j], format = "d", width = col_widths[6]),
            formatC(results$Nh1[j], format = "d", width = col_widths[7])
          )
        } else {
          row_vals <- c(
            formatC(j, width = col_widths[1], format = "d"),
            formatC(results$h0[j], format = "f", digits = 3, width = col_widths[2]),
            formatC(results$h1[j], format = "f", digits = 3, width = col_widths[3]),
            formatC(results$Nh0[j], format = "d", width = col_widths[4]),
            formatC(results$Nh1[j], format = "d", width = col_widths[5])
          )
        }
        cat(paste(row_vals, collapse = ""), "\n")
      }
    }
    cat(strrep("=", sum(col_widths)), "\n")
  }
}

