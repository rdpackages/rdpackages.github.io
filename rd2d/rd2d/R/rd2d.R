################################################################################
#' @title Two-Dimensional Local Polynomial Regression Discontinuity Design
#'
#' @description
#' \code{rd2d} implements bivariate local polynomial boundary regression discontinuity (RD) point estimators with robust bias-corrected pointwise confidence intervals and 
#' uniform confidence bands, developed in Cattaneo, Titiunik and Yu (2025). For robust bias-correction, see Calonico, Cattaneo and Titiunik (2014).
#' 
#' Companion commands are: \code{rdbw2d} for data-driven bandwidth selection.
#' 
#' For other packages of RD designs, visit
#' <https://rdpackages.github.io/>
#'
#' @param Y Dependent variable; a numeric vector of length \eqn{N}, where \eqn{N} is the sample size.
#' @param X Bivariate running variable (a.k.a score variable); a numeric matrix or data frame of dimension \eqn{N \times 2}, with each row \eqn{\mathbf{X}_i = (X_{1i}, X_{2i})}.
#' @param t Treatment indicator; a logical or binary vector indicating treatment assignment (\eqn{t_i = 1} if treated, \eqn{t_i = 0} otherwise).
#' @param b Evaluation points; a matrix or data frame specifying boundary points \eqn{\mathbf{b}_j = (b_{1j}, b_{2j})}, of dimension \eqn{J \times 2}.
#' @param h Bandwidths. Either a positive scalar (same bandwidth for all dimensions and groups), or a matrix/data frame of size \eqn{J \times 4}, corresponding to \eqn{h_{\text{control},1}}, 
#' \eqn{h_{\text{control},2}}, \eqn{h_{\text{treated},1}}, \eqn{h_{\text{treated},2}} at each evaluation point. If not specified, bandwidth \code{h} is computed by the companion command 
#' \code{rdbw2d}. Default is \code{h = NULL}.
#' @param deriv The order of the derivatives of the regression functions to be estimated; a numeric vector of length 2 specifying the number of derivatives in each coordinate (e.g., \eqn{c(1,2)} corresponds to \eqn{\partial_1 \partial_2^2}).
#' @param tangvec Tangent vectors; a matrix or data frame of dimension \eqn{J \times 2} specifying directional derivatives. Overrides \code{deriv} if provided.
#' @param p Polynomial order for point estimation (\eqn{p = 1} by default).
#' @param q Polynomial order for robust confidence interval construction. Must satisfy \eqn{q \geq p}; default is \eqn{q = p + 1}.
#' @param kernel Kernel function to use. Options are \code{"unif"}, \code{"uniform"} (uniform), \code{"triag"}, \code{"triangular"} (triangular, default), and \code{"epan"}, \code{"epanechnikov"} (Epanechnikov).
#' @param kernel_type Kernel structure. Either \code{"prod"} for product kernels or \code{"rad"} for radial kernels.
#' @param vce Variance-covariance estimation method. Options are:
#' \itemize{
#' \item \code{"hc0"}: heteroskedasticity-robust plug-in residual variance estimator without small-sample adjustment.
#' \item \code{"hc1"}: heteroskedasticity-robust plug-in residual variance estimator with HC1 small-sample adjustment (default).
#' \item \code{"hc2"}: heteroskedasticity-robust plug-in residual variance estimator with HC2 adjustment.
#' \item \code{"hc3"}: heteroskedasticity-robust plug-in residual variance estimator with HC3 adjustment.
#' }
#' Default is \code{"hc1"}.
#' @param masspoints Handling of mass points in the running variable. Options are:
#' \itemize{
#' \item \code{"check"}: detects presence of mass points and reports the number of unique observations (default).
#' \item \code{"adjust"}: adjusts preliminary bandwidths to ensure a minimum number of unique observations within each side of the cutoff.
#' \item \code{"off"}: ignores presence of mass points.
#' }
#' @param C Cluster ID variable used for cluster-robust variance estimation with degrees-of-freedom weights.Default is \code{C = NULL}.
#' @param level Nominal confidence level for intervals/bands, between 0 and 100 (default is 95).
#' @param cbands Logical. If \code{TRUE}, also compute uniform confidence bands (default is \code{FALSE}).
#' @param side Type of confidence interval. Options: \code{"two"} (two-sided, default), \code{"left"} (left tail), or \code{"right"} (right tail).
#' @param repp Number of repetitions for critical value simulation (used in uniform confidence bands). Default is 1000.
#' @param bwselect Bandwidth selection strategy. Options: 
#' \itemize{
#' \item \code{"mserd"}. One common MSE-optimal bandwidth selector for the boundary RD treatment effect estimator for each evaluation point (default). 
#' \item \code{"imserd"}. IMSE-optimal bandwidth selector for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"msetwo"}. Two different MSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator for each evaluation point. 
#' \item \code{"imsetwo"}. Two IMSE-optimal bandwidth selectors (control and treatment) for the boundary RD treatment effect estimator based on all evaluation points.
#' \item \code{"user provided"}. User-provided bandwidths. If \code{h} is not \code{NULL}, then \code{bwselect} is overwritten to \code{"user provided"}.
#' }
#' @param method Bandwidth selection method for bias estimator based on local polynomials. Either \code{"dpi"} (default) for data-driven plug-in MSE optimal bandwidth selector or \code{"rot"} for rule-of-thumb bandwidth selector.
#' @param bwcheck If a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least \code{bwcheck} unique observations are used. Default is \code{50 + p + 1}.
#' @param scaleregul Scaling factor for the regularization term in bandwidth selection. Default is 3.
#' @param scalebiascrct Scaling factor used for bias correction based on higher order expansions. Default is 1.
#' @param stdvars Logical. If TRUE, the running variables \eqn{X_{1i}} and \eqn{X_{2i}} are standardized before computing the bandwidths. Default is \code{FALSE}. Standardization only affects automatic bandwidth selection if bandwidths are not manually provided via \code{h}.
#'
#' @return An object of class \code{"rd2d"}, a list with components:
#' \describe{
#'   \item{\code{main}}{A data frame with point estimates, variances, p-values, confidence intervals, confidence bands, and bandwidths at each evaluation point. 
#'     \describe{
#'       \item{\code{b1}, \code{b2}}{First and second coordinate of evaluation points \eqn{\mathbf{b} = (b_1,b_2)}.}
#'       \item{\code{Est.p}}{Point estimate \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'       \item{\code{Var.p}}{Variance of estimate \eqn{\widehat{\tau}_p(\mathbf{b})}.}
#'       \item{\code{Est.q}}{Bias-corrected point estimate \eqn{\widehat{\tau}_q(\mathbf{b})}.}
#'       \item{\code{Var.q}}{Variance of bias-corrected estimate \eqn{\widehat{\tau}_q(\mathbf{b})}.}
#'       \item{\code{p-value}}{P-value based on t-statistic with bias-corrected estimate.}
#'       \item{\code{CI.lower}, \code{CI.upper}}{Pointwise confidence intervals.}
#'       \item{\code{CB.lower}, \code{CB.upper}}{Uniform confidence bands if computed.}
#'       \item{\code{h01}, \code{h02}, \code{h11}, \code{h12}}{Bandwidths used in each coordinate and group.}
#'       \item{\code{Nh0}, \code{Nh1}}{Effective sample size on each side of the cutoff.}
#'     }
#'   }
#'   \item{\code{main.A0}}{Same structure as \code{main} but for control group outcomes.}
#'   \item{\code{main.A1}}{Same structure as \code{main} but for treated group outcomes.}
#'   \item{\code{tau.hat}}{Estimated treatment effect at each evaluation point.}
#'   \item{\code{se.hat}}{Standard errors corresponding to estimates.}
#'   \item{\code{cov.us}}{Covariance matrix used for uniform bands.}
#'   \item{\code{cb}}{List with critical values, pointwise, and uniform intervals.}
#'   \item{\code{pvalues}}{Two-sided p-values based on bias-corrected estimates.}
#'   \item{\code{opt}}{List of options used in the function call.}
#' }
#'   
#' @seealso \code{\link{rdbw2d}}, \code{\link{print.rd2d}}, \code{\link{summary.rd2d}}
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @references
#' \itemize{
#' \item{\href{https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}{Calonico, S., M. D. Cattaneo, and R. Titiunik. (2014).}
#' Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs}
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
#' # Estimate treatment effect using rd2d
#' result <- rd2d(Y, X, t, b, cbands = TRUE)
#' print(result)
#' summary(result)
#' @export
rd2d <- function(Y, X, t, b, h = NULL, deriv = c(0,0), tangvec = NULL,
                 p = 1, q = 2, kernel = c("tri","triangular","epa","epanechnikov","uni","uniform","gau","gaussian"),
                 kernel_type = c("prod","rad"), vce = c("hc1","hc0","hc2","hc3"),
                 masspoints = c("check", "adjust", "off"),C = NULL,
                 level = 95, cbands = TRUE, side = c("two", "left", "right"), repp = 1000,
                 bwselect = c("mserd", "imserd", "msetwo", "imsetwo", "user provided"),
                 method = c("dpi", "rot"), bwcheck = 50 + p + 1,
                 scaleregul = 3, scalebiascrct = 1, stdvars = TRUE){

  ######################## Input error handling ################################

  kernel <- match.arg(kernel)
  kernel_type <- match.arg(kernel_type)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  side <- match.arg(side)
  bwselect <- match.arg(bwselect)
  method <- match.arg(method)
  
  d <- t # renaming the variable

  exit <- 0

  # Check Y, X, d lengths
  if (length(Y) != length(d) || length(Y) != nrow(X)) {
    print("Y, d, and rows of X must have the same length")
    exit <- 1
  }

  # X must have 2 columns
  if (ncol(X) != 2) {
    print("X must have exactly 2 columns")
    exit <- 1
  }

  # d must be logical or contain only 0 and 1
  if (!(is.logical(d) || all(d %in% c(0, 1)))) {
    print("d must be a logical vector or a numeric vector containing only 0 and 1")
    exit <- 1
  }

  # b must be a matrix with 2 columns
  if (!(is.matrix(b) || is.data.frame(b)) || ncol(b) != 2) {
    print("b must be a matrix with 2 columns")
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
               nrow(h) != nrow(b) || ncol(h) != 4) {
      print("If h is not a scalar, it must be a matrix or data frame with the same number of rows as b and 4 columns")
      exit <- 1
    }
  }

  # deriv must be a numeric vector of length 2, and deriv[1] + deriv[2] <= p
  if (!is.numeric(deriv) || length(deriv) != 2) {
    print("deriv must be a numeric vector of length 2")
    exit <- 1
  } else if (sum(deriv) > p) {
    print("Sum of deriv components must be less than or equal to polynomial order p")
    exit <- 1
  }

  # tangvec, if provided, must be matrix/data.frame with same nrow as b and 2 columns
  if (!is.null(tangvec)) {
    if (!(is.matrix(tangvec) || is.data.frame(tangvec)) ||
        nrow(tangvec) != nrow(b) || ncol(tangvec) != 2) {
      warning("tangvec must be a matrix or data frame with same number of rows as b and 2 columns")
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

  if (q < p){
    print("Parameter q must be no smaller than p. Please provide valid inputs.")
    exit <- 1
  }

  if (!is.null(C) && !(vce %in% c("hc0", "hc1"))) {
    warning("When C is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }

  if (exit>0) stop()

  ############################ Data preparation ################################

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

  # Check for mass points

  M <- NULL; M.0 <- NULL; M.1 <- NULL


  # Only check for mass points if inputting bivariate coordinates.
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

  ################################ Bandwidth ###################################

  if (is.null(h)){
    # if (bwselect == "mserd" | bwselect == "imserd"){
    #   bws <- (rdbw2d(dat, eval, o, p, deriv, kernel, bwselect, method, vce, bwcheck, masspoints, C, scaleregul, bydist, inputdist))$bws
    #   hgrid <- bws$hMSE.0
    #   hgrid.1 <- NULL
    # } else if (bwselect == "msetwo" | bwselect == "imsetwo"){
    #   bws <- (rdbw2d(dat, eval, o, p, deriv, kernel, bwselect, method, vce, bwcheck, masspoints, C, scaleregul, bydist, inputdist))$bws
    #   hgrid <- bws$hMSE.0
    #   hgrid.1 <- bws$hMSE.1
    # }
    bws <- rdbw2d(Y = Y, X = X, t = d, b = b, p = p, deriv = deriv, tangvec = tangvec,
                  kernel = kernel, kernel_type = kernel_type,
                  bwselect = bwselect, method = method, vce = vce,
                  bwcheck = bwcheck, masspoints = masspoints,
                  C = C, scaleregul = scaleregul, scalebiascrct = scalebiascrct,
                  stdvars = stdvars)
    bws <- bws$bws
    hgrid <- cbind(bws[,3],bws[,4])
    hgrid.1 <- cbind(bws[,5],bws[,6])
  } else {
    bwselect <- "user provided"
    # standardize bandwidth
    if (length(h) == 1){
      hgrid <- matrix(h, nrow = neval, ncol = 2)
      hgrid.1 <- matrix(h, nrow = neval, ncol = 2)
    } else {
      hgrid <- cbind(h[,1],h[,2])
      hgrid.1 <- cbind(h[,3],h[,4])
    }
  }

  ###################### Point estimation and inference ########################

  count.q <- factorial(q + 2)/(factorial(q) * 2)
  count.p <- factorial(p + 2)/(factorial(p) * 2)

  e_deriv.q <- matrix(0, nrow = neval, ncol = count.q)
  e_deriv.q[,c(1:count.p)] <- e_deriv

  rdfit.p <- rd2d_fit_v2(dat, eval, e_deriv, deriv, p, hgrid, hgrid.1, kernel, kernel_type, vce, masspoints, C, bwcheck, unique)
  tau.hat.p <- rdfit.p$mu.1 - rdfit.p$mu.0
  se.hat.p <- sqrt(rdfit.p$se.0^2 + rdfit.p$se.1^2)
  h.01.p <- rdfit.p$h.0.x
  h.02.p <- rdfit.p$h.0.y
  h.11.p <- rdfit.p$h.1.x
  h.12.p <- rdfit.p$h.1.y
  eN.0.p <- rdfit.p$eN.0
  eN.1.p <- rdfit.p$eN.1

  rdfit.q <- rd2d_fit_v2(dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel, kernel_type, vce, masspoints, C, bwcheck, unique)
  tau.hat.q <- rdfit.q$mu.1 - rdfit.q$mu.0
  se.hat.q <- sqrt(rdfit.q$se.0^2 + rdfit.q$se.1^2)
  h.01.q <- rdfit.q$h.0.x
  h.02.q <- rdfit.q$h.0.y
  h.11.q <- rdfit.q$h.1.x
  h.12.q <- rdfit.q$h.1.y
  eN.0.q <- rdfit.q$eN.0
  eN.1.q <- rdfit.q$eN.1

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

  # Covariance

  cov.hat.q <- NA
  cb.hat.q <- list(CI.l = CI.lower, CI.r = CI.upper, CB.l = rep(NA, length(CI.lower)), CB.r = rep(NA, length(CI.lower)))
  CB.lower <- NA
  CB.upper <- NA
  if (cbands){
    cov.hat.q <- rdbw2d_cov(dat, eval, e_deriv.q, deriv, q, hgrid, hgrid.1, kernel, kernel_type, vce, C)
    cb.hat.q <- rd2d_cb(tau.hat.q, cov.hat.q, repp, side, level)
    CB.lower <- cb.hat.q$CB.l
    CB.upper <- cb.hat.q$CB.r
  }

  clustered <- !is.null(C)

  ################################## Output ####################################

  main <- cbind(b[,1], b[,2], tau.hat.p, se.hat.p, tau.hat.q, se.hat.q, zvalues, pvalues,
                CI.lower, CI.upper, CB.lower, CB.upper, hgrid[,1], hgrid[,2],
                hgrid.1[,1], hgrid.1[,2], eN.0.p, eN.1.p)
  main <- as.data.frame(main)
  colnames(main) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q", "z", "P>|z|",
                      "CI.lower","CI.upper","CB.lower", "CB.upper", "h01", "h02",
                      "h11", "h12", "Nh0", "Nh1")

  main.A0 <- cbind(b[,1], b[,2], rdfit.p$mu.0, rdfit.p$se.0, rdfit.q$mu.0,rdfit.q$se.0, hgrid[,1], hgrid[,2], eN.0.p)
  main.A0 <- as.data.frame(main.A0)
  colnames(main.A0) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q","h01", "h02","Nh0")

  main.A1 <- cbind(b[,1], b[,2], rdfit.p$mu.1, rdfit.p$se.1, rdfit.q$mu.1,rdfit.q$se.1, hgrid.1[,1], hgrid.1[,2], eN.1.p)
  main.A1 <- as.data.frame(main.A1)
  colnames(main.A1) <- c("b1","b2","Est.p","Se.p","Est.q","Se.q","h11", "h12","Nh1")

  rdmodel <- "rd2d"

  out <- list(main = main, main.A0 = main.A0, main.A1 = main.A1,
              opt=list(b = b, deriv = deriv, tangvec = tangvec, p = p, q = q, kernel=kernel.type, kernel_type = kernel_type, N=N, N.0 = N.0,
                       N.1 = N.1, M = M, M.0 = M.0, M.1 = M.1, neval=neval, bwselect = bwselect, method = method,
                       vce = vce, bwcheck = bwcheck, masspoints = masspoints, C = C, clustered = clustered,
                       scaleregul = scaleregul, scalebiascrct = scalebiascrct, stdvars = stdvars,
                       level = level, repp = repp, side = side,cbands = cbands,
                       h01 = hgrid[,1], h02 = hgrid[,2], h11 = hgrid.1[,1], h12 = hgrid.1[,2],
                       Nh0 = eN.0.p, Nh1 = eN.1.p),
              tau.hat = tau.hat.p, tau.hat.q = tau.hat.q, se.hat = se.hat.p, se.hat.q = se.hat.q, cov.us=cov.hat.q, cb = cb.hat.q, pvalues = pvalues,
              zvalues = zvalues,rdmodel = rdmodel)
  out$call   <- match.call()
  class(out) <- "rd2d"

  return(out)
}

################################################################################
#' Print Method for 2D Local Polynomial RD Estimation
#' @description
#' Prints the results of a 2D local polynomial regression discontinuity (RD) estimation, as obtained from \code{\link{rd2d}}.
#'
#' @param x An object of class \code{rd2d}, returned by \code{\link{rd2d}}.
#' @param ... Additional arguments passed to the method (currently ignored).
#'
#' @return
#' No return value. This function is called for its side effects, which are to print the \code{\link{rd2d}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso
#' \code{\link{rd2d}} for conducting 2D local polynomial RD estimation.
#'
#' Supported methods: \code{\link{print.rd2d}}, \code{\link{summary.rd2d}}.
#'
#' @export

print.rd2d <- function(x,...) {

  cat(paste(x$rdmodel, "\n", sep = ""))
  cat(paste("\n", sep = ""))

  # Format and print the vector as "(x, y)"
  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  cat(sprintf("BW type.               %s\n", paste(paste(x$opt$bwselect, x$opt$method, sep = "-"), ifelse(x$opt$stdvar, "-std", ""),sep = "")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), x$opt$kernel_type, sep = "-")))
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  # cat(sprintf("Standardization        %s\n", ifelse(x$opt$stdvars, "on", "off")))
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Estimand (deriv)       %-10d   %-10d\n", x$opt$deriv[1], x$opt$deriv[2]))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
  if (x$opt$masspoints == "check" | x$opt$masspoints == "adjust") {
    cat(sprintf("Unique Obs.            %-10d   %-10d\n", x$opt$M.0, x$opt$M.1))
  }
  cat("\n")

}


################################################################################
#' Summary Method for 2D Local Polynomial RD Estimation
#'
#' @description
#' Summarizes estimation and bandwidth results from a 2D local polynomial regression discontinuity (RD) design, as produced by \code{\link{rd2d}}.
#'
#' @param object An object of class \code{rd2d}, typically returned by \code{\link{rd2d}}.
#' @param ... Optional arguments. Supported options include:
#'   \itemize{
#'     \item \code{CBuniform}: Logical. If \code{TRUE}, displays uniform confidence bands;
#'       if \code{FALSE} (default), displays pointwise confidence intervals.
#'     \item \code{subset}: Integer vector of indices of evaluation points to display.
#'       Defaults to all evaluation points.
#'     \item \code{output}: Character. Use \code{"main"} to display estimation results,
#'       or \code{"bw"} to display bandwidth information. Default is \code{"main"}.
#'     \item \code{sep}: Integer vector of length three. Controls spacing in the output.
#'       \code{sep[1]} controls spacing for the columns of bandwidths, estimation,
#'       z-value, and p-value in the \code{"main"} table.
#'       \code{sep[2]} controls spacing for the confidence interval (confidence bands)
#'       in the \code{"main"} table.
#'       \code{sep[3]} controls spacing for the columns in the \code{"bw"} table.
#'       Default is \code{c(7, 17, 8)}.
#'   }
#'
#' @return No return value. This function is called for its side effects: it prints a formatted summary of \code{\link{rd2d}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso \code{\link{rd2d}} for estimation using 2D local polynomial RD design.
#'
#' Supported methods: \code{\link{print.rd2d}}, \code{\link{summary.rd2d}}.
#'
#' @export

summary.rd2d <- function(object, ...) {

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

  # Format and print the vector as "(x, y)"
  cat(sprintf("Number of Obs.         %d\n", x$opt$N))
  cat(sprintf("BW type.               %s\n", paste(paste(x$opt$bwselect, x$opt$method, sep = "-"), ifelse(x$opt$stdvar, "-std", ""),sep = "")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), x$opt$kernel_type, sep = "-")))
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  # cat(sprintf("Standardization        %s\n", ifelse(x$opt$stdvars, "on", "off")))
  cat("\n")
  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Estimand (deriv)       %-10d   %-10d\n", x$opt$deriv[1], x$opt$deriv[2]))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat(sprintf("Order rbc. (q)         %-10d   %-10d\n", x$opt$q, x$opt$q))
  if (x$opt$masspoints == "check" | x$opt$masspoints == "adjust") {
    cat(sprintf("Unique Obs.            %-10d   %-10d\n", x$opt$M.0, x$opt$M.1))
  }
  cat("\n")

  # Prepare data for the bottom table

  results <- data.frame(
    b1 = x$opt$b[, 1],
    b2 = x$opt$b[, 2],
    Coef = x$tau.hat,
    Zvalues = x$zvalues,
    Pvalues = x$pvalues,
    CILower = if (CBuniform) x$cb$CB.l else x$cb$CI.l,
    CIUpper = if (CBuniform) x$cb$CB.r else x$cb$CI.r,
    h01 = x$opt$h01,
    h02 = x$opt$h02,
    h11 = x$opt$h11,
    h12 = x$opt$h12,
    Nh0 = x$opt$Nh0,
    Nh1 = x$opt$Nh1
  )

  if (output == "main"){
    headers <- c("ID", "b1", "b2", "Est.", "z", "P>|z|", sprintf("%d%% CI", x$opt$level))
    if (CBuniform){
      headers[length(headers)] <- sprintf("%d%% Unif. CB", x$opt$level)
    }

    col_widths <- c(4, sep[1], sep[1], sep[1], sep[1], sep[1], sep[2])

    # Format and print header row
    cat(paste(rep("=", sum(col_widths) + 2 * (length(headers) - 1)), collapse = ""), "\n")
    formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
    cat(paste(formatted_headers, collapse = "  "), "\n")
    cat(paste(rep("=", sum(col_widths) + 2 * (length(headers) - 1)), collapse = ""), "\n")

    neval <- nrow(results)
    if (is.null(subset)){
      subset <- seq_len(neval)
    } else{
      # input error handling
      if (!all(subset %in% seq_len(neval))) {
        warning("Invalid subset provided. Resetting to default: 1:neval")
        subset <- seq_len(neval)
      }
    }

    # Print each row of the results
    for (i in 1:nrow(results)) {

      index_formatted <- formatC(i, width = col_widths[1], format = "d")

      b1_formatted <- ifelse(
        is.na(results$b1[i]),
        formatC("NA", width = col_widths[2], format = "s"),
        formatC(results$b1[i], format = "f", digits = 3, width = col_widths[2])
      )

      b2_formatted <- ifelse(
        is.na(results$b2[i]),
        formatC("NA", width = col_widths[3], format = "s"),
        formatC(results$b2[i], format = "f", digits = 3, width = col_widths[3])
      )

      coef_formatted <- ifelse(
        is.na(results$Coef[i]),
        formatC("NA", width = col_widths[4], format = "s"),
        formatC(results$Coef[i], format = "f", digits = 4, width = col_widths[4])
      )

      zvalues_formatted <- ifelse(
        is.na(results$Zvalues[i]),
        formatC("NA", width = col_widths[5], format = "s"),
        formatC(results$Zvalues[i], format = "f", digits = 4, width = col_widths[5])
      )

      pvalues_formatted <- ifelse(
        is.na(results$Pvalues[i]),
        formatC("NA", width = col_widths[6], format = "s"),
        formatC(results$Pvalues[i], format = "f", digits = 4, width = col_widths[6])
      )

      ci_formatted <- ifelse(
        is.na(results$CILower[i]) | is.na(results$CIUpper[i]),
        formatC("NA", width = col_widths[7], format = "s"),
        formatC(paste0("[", formatC(results$CILower[i], format = "f", digits = 4),
                       ", ", formatC(results$CIUpper[i], format = "f", digits = 4), "]"),
                width = col_widths[7], format = "s")
      )

      # Print
      row_vals <- c(index_formatted,
                    b1_formatted, b2_formatted, coef_formatted, zvalues_formatted,pvalues_formatted,
                    ci_formatted
      )
      if (i %in% subset) cat(paste(row_vals, collapse = "  "), "\n")
    }
    cat(paste(rep("=", sum(col_widths) + 2 * (length(headers) - 1)), collapse = ""), "\n")
  } else {
    headers <- c("ID", "b1", "b2", "h01", "h02", "h11", "h12", "Nh0", "Nh1")
    col_widths <- c(4, sep[3], sep[3], sep[3], sep[3], sep[3], sep[3], sep[3], sep[3])

    cat(strrep("=", sum(col_widths)), "\n")
    group_headers <- c(
      formatC("        Bdy Points", width = col_widths[1] + col_widths[2] + col_widths[3], format = "s", flag = "-"),
      formatC("      BW Control", width = col_widths[4] + col_widths[5], format = "s", flag = "-"),
      formatC("    BW Treatment", width = col_widths[6] + col_widths[7], format = "s", flag = "-"),
      formatC("       Eff. N", width = col_widths[8] + col_widths[9], format = "s", flag = "-")
    )
    cat(paste(group_headers, collapse = ""), "\n")

    # Format and print header row
    formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
    cat(paste(formatted_headers, collapse = ""), "\n")
    cat(strrep("=", sum(col_widths)), "\n")

    neval <- nrow(results)
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
    for (j in 1:nrow(results)) {
      index <- formatC(j, width = col_widths[1], format = "d")
      bdy1 <- ifelse(is.na(results[j, "b1"]),
                     formatC("NA", width = col_widths[2], format = "s"),
                     formatC(results[j, "b1"], format = "f", digits = 3, width = col_widths[2]))
      bdy2 <- ifelse(is.na(results[j, "b2"]),
                     formatC("NA", width = col_widths[3], format = "s"),
                     formatC(results[j, "b2"], format = "f", digits = 3, width = col_widths[3]))
      control1 <- ifelse(is.na(results[j, "h01"]),
                         formatC("NA", width = col_widths[4], format = "s"),
                         formatC(results[j, "h01"], format = "f", digits = 3, width = col_widths[4]))
      control2 <- ifelse(is.na(results[j, "h02"]),
                         formatC("NA", width = col_widths[5], format = "s"),
                         formatC(results[j, "h02"], format = "f", digits = 3, width = col_widths[5]))
      treatment1 <- ifelse(is.na(results[j, "h11"]),
                           formatC("NA", width = col_widths[6], format = "s"),
                           formatC(results[j, "h11"], format = "f", digits = 3, width = col_widths[6]))
      treatment2 <- ifelse(is.na(results[j, "h12"]),
                           formatC("NA", width = col_widths[7], format = "s"),
                           formatC(results[j, "h12"], format = "f", digits = 3, width = col_widths[7]))
      Nh0 <- ifelse(is.na(results[j, "Nh0"]),
                           formatC("NA", width = col_widths[8], format = "s"),
                           formatC(results[j, "Nh0"], format = "d", width = col_widths[8]))
      Nh1 <- ifelse(is.na(results[j, "Nh1"]),
                           formatC("NA", width = col_widths[9], format = "s"),
                           formatC(results[j, "Nh1"], format = "d", width = col_widths[9]))
      # Combine formatted values and print the row
      row_vals <- c(index, bdy1, bdy2, control1, control2, treatment1, treatment2, Nh0, Nh1)
      if (j %in% subset) cat(paste(row_vals, collapse = ""), "\n")
    }
    cat(strrep("=", sum(col_widths)), "\n")
  }

}



