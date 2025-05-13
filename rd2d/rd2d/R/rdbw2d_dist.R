# This file is for bandwidth selection for local polynomial estimator based on 2d location data.
#' @title Bandwidth Selection for Distance-Based RD Designs
#' @description
#' \code{rdbw2d.dist} implements bandwidth selector for distance-based local polynomial boundary regression discontinuity (RD) point estimators with robust bias-corrected pointwise confidence intervals and 
#' uniform confidence bands, developed in Cattaneo, Titiunik and Yu (2025). For robust bias-correction, see Calonico, Cattaneo and Titiunik (2014).
#'
#' @param Y Dependent variable; a numeric vector of length \eqn{N}, where \eqn{N} is the sample size.
#' @param D Distance-based scores \eqn{\mathbf{D}_i=(\mathbf{D}_{i}(\mathbf{b}_1),\cdots,\mathbf{D}_{i}(\mathbf{b}_J))}; dimension is \eqn{N \times J} where \eqn{N} = sample size and \eqn{J} = number of cutoffs; 
#' non-negative values means data point in treatment group and negative values means data point in control group.
#' @param b Optional evaluation points; a matrix or data frame specifying boundary points \eqn{\mathbf{b}_j = (b_{1j}, b_{2j})}, dimension \eqn{J \times 2}.
#' @param p Polynomial order for point estimation. Default is \code{p = 1}.
#' @param kink Logical; whether to apply kink adjustment. Options: \code{"on"} (default) or \code{"off"}.
#' @param kernel Kernel function to use. Options are \code{"unif"}, \code{"uniform"} (uniform), \code{"triag"}, \code{"triangular"} (triangular, default), and \code{"epan"}, \code{"epanechnikov"} (Epanechnikov).
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
#' @return An object of class \code{"rdbw2d.dist"}, containing:
#' \describe{
#'   \item{\code{bws}}{Data frame of optimal bandwidths for each evaluation point:
#'     \describe{
#'       \item{\code{b1}}{First coordinate of the evaluation point \eqn{b1}.}
#'       \item{\code{b2}}{Second coordinate of the evaluation point \eqn{b2}.}
#'       \item{\code{h0}}{Bandwidth for observations with distance \eqn{D_{i}(\mathbf{b}) < 0}.}
#'       \item{\code{h1}}{Bandwidth for observations with distance \eqn{D_{i}(\mathbf{b}) \geq 0}.}
#'       \item{\code{Nh0}}{Effective sample size for \eqn{D_{i}(\mathbf{b}) < 0}.}
#'       \item{\code{Nh1}}{Effective sample size for \eqn{D_{i}(\mathbf{b}) \geq 0}.}
#'     }
#'   }
#'   \item{\code{mseconsts}}{Data frame of intermediate bias and variance constants used for MSE/IMSE calculations.}
#'   \item{\code{opt}}{A list of options and settings used in estimation, including \code{p}, \code{kernel}, sample size \eqn{N}, and user-specified choices.}
#' }
#'
#' @seealso \code{\link{rd2d.dist}}, \code{\link{rd2d}}, \code{\link{summary.rdbw2d.dist}}, \code{\link{print.rdbw2d.dist}}
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
#' bws <- rdbw2d.dist(y, D, b = eval)
#'
#' # View the estimation results
#' print(bws)
#' summary(bws)
#' @export



rdbw2d.dist <- function(Y, D, b = NULL, p = 1, kink = c("off", "on"),
                   kernel = c("tri","triangular", "epa","epanechnikov","uni","uniform","gau","gaussian"),
                   bwselect = c("mserd", "imserd", "msetwo", "imsetwo"),
                   vce = c("hc1","hc0","hc2","hc3"),
                   bwcheck = 20 + p + 1, masspoints = c("check","adjust","off"),
                   C = NULL, scaleregul = 1, cqt = 0.5){

  # Input error handling

  bwselect <- match.arg(bwselect)
  kernel <- match.arg(kernel)
  vce <- match.arg(vce)
  masspoints <- match.arg(masspoints)
  kink <- match.arg(kink)
  rot <- NULL

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

  if (!is.null(C) && !(vce %in% c("hc0", "hc1"))) {
    warning("When C is specified, vce must be 'hc0' or 'hc1'. Resetting vce to 'hc1'.")
    vce <- "hc1"
  }

  if (exit>0) stop()

  # Data Cleaning

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

  e_deriv <- matrix(0, nrow = neval, ncol = p+1)
  e_deriv[,1] <- 1

  kernel.type <- "Epanechnikov"
  if (kernel=="triangular"   | kernel=="tri") kernel.type <- "Triangular"
  if (kernel=="uniform"      | kernel=="uni") kernel.type <- "Uniform"
  if (kernel=="gaussian"     | kernel=="gau") kernel.type <- "Gaussian"

  # Store variance and bias constants for IMSE

  bconst <- rep(NA, neval)
  vconst <- rep(NA, neval)

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

  results <- rdbw2d_dist_bw(Y = Y, D = D, p = p, kernel = kernel, deriv.vec = e_deriv, rot = rot,
                            vce = vce, C = C, bwcheck = bwcheck, scaleregul = scaleregul, cqt = cqt)

  if (bwselect == "mserd"){
    deriv.sum <- 0
    hn.grid <- ( (2 + 2 * deriv.sum) * (results$v.0   + results$v.1) /
              ( (2 * p + 2 - 2 * deriv.sum) * ( (results$b.0 - results$b.1)^2 +
              scaleregul * results$r.0 + scaleregul * results$r.1) ) )^(1/(2 * p + 4))
    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid <- pmax(hn.grid, results$bw.min.0, results$bw.min.1)
      hn.grid <- pmin(hn.grid, results$bw.max.0, results$bw.max.1)
    }
    results$h.0 <- hn.grid
    results$h.1 <- hn.grid
  }

  if (bwselect == "msetwo"){
    deriv.sum <- 0
    hn.grid.0 <- ( (2 + 2 * deriv.sum) * results$v.0 /
                   ( (2 * p + 2 - 2 * deriv.sum) * ( results$b.0^2 + scaleregul * results$r.0) ) )^(1/(2 * p + 4))
    hn.grid.1 <- ( (2 + 2 * deriv.sum) * results$v.1 /
                   ( (2 * p + 2 - 2 * deriv.sum) * ( results$b.1^2 + scaleregul * results$r.1) ) )^(1/(2 * p + 4))
    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid.0 <- pmax(hn.grid.0, results$bw.min.0)
      hn.grid.0 <- pmin(hn.grid.0, results$bw.max.0)
      hn.grid.1 <- pmax(hn.grid.1, results$bw.min.1)
      hn.grid.1 <- pmin(hn.grid.1, results$bw.max.1)
    }
    results$h.0 <- hn.grid.0
    results$h.1 <- hn.grid.1
  }

  if (bwselect == "imserd"){
    V.V <- mean(results$v.0) + mean(results$v.1)
    B.B <- mean( (results$b.0 - results$b.1)^2 + scaleregul * results$r.0 + scaleregul * results$r.1)
    deriv.sum <- 0
    hIMSE <- ((2 + 2 * deriv.sum) * V.V / ( (2 * p + 2)* B.B ) )^(1/(2 * p + 4))
    hn.grid <- rep(hIMSE, nrow(results))
    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid <- pmax(hn.grid, results$bw.min.0, results$bw.min.1)
      hn.grid <- pmin(hn.grid, results$bw.max.0, results$bw.max.1)
    }
    results$h.0 <- hn.grid
    results$h.1 <- hn.grid
  }

  if (bwselect == "imsetwo"){
    V.V.0 <- mean(results$v.0)
    V.V.1 <- mean(results$v.1)
    B.B.0 <- mean( results$b.0^2 + scaleregul * results$r.0)
    B.B.1 <- mean( results$b.1^2 + scaleregul * results$r.1)
    deriv.sum <- 0
    hIMSE.0 <- ((2 + 2 * deriv.sum) * V.V.0 / ( (2 * p + 2)* B.B.0 ) )^(1/(2 * p + 4))
    hIMSE.1 <- ((2 + 2 * deriv.sum) * V.V.1 / ( (2 * p + 2)* B.B.1 ) )^(1/(2 * p + 4))
    hn.grid.0 <- rep(hIMSE.0, nrow(results))
    hn.grid.1 <- rep(hIMSE.1, nrow(results))

    if (!is.null(bwcheck)) { # Bandwidth restrictions
      hn.grid.0 <- pmax(hn.grid.0, results$bw.min.0)
      hn.grid.0 <- pmin(hn.grid.0, results$bw.max.0)
      hn.grid.1 <- pmax(hn.grid.1, results$bw.min.1)
      hn.grid.1 <- pmin(hn.grid.1, results$bw.max.1)
    }
    results$h.0 <- hn.grid.0
    results$h.1 <- hn.grid.1
  }

  # kink adjustment

  if (kink =="on"){
    if (bwselect == "mserd" | bwselect == "imserd"){
      results$h.0 <- results$h.0 * M.vec^(-1/4) / M.vec^(-1/(2 * p + 4))
      results$h.1 <- results$h.1 * M.vec^(-1/4) / M.vec^(-1/(2 * p + 4))
    }
    if (bwselect == "msetwo" | bwselect == "imsetwo"){
      results$h.0 <- results$h.0 * M.0.vec^(-1/4) / M.0.vec^(-1/(2 * p + 4))
      results$h.1 <- results$h.1 * M.1.vec^(-1/4) / M.1.vec^(-1/(2 * p + 4))
    }
  }

  # Outputs

  bws <- results[,c("h.0", "h.1")]
  bws <- cbind(eval, bws)
  colnames(bws) <- c("b1","b2","h0", "h1")

  out        <- list(bws = bws, mseconsts = results,
                     opt = list(N=N, N.0 = N.0, N.1 = N.1, M = M.vec, M.0 = M.0.vec,
                                M.1 = M.1.vec, neval=neval, p=p, b = eval,
                                kernel=kernel.type, kink = kink,
                                bwselect=bwselect, bwcheck = bwcheck,
                                C= C,
                                vce = vce, masspoints = masspoints,
                                scaleregul = scaleregul, cqt = cqt))
  out$call   <- match.call()
  class(out) <- "rdbw2d.dist"
  return(out)
}

# Print method

################################################################################
#' Print Method for Bandwidth Selection (Distance-Based) in 2D Local Polynomial RD Design
#'
#' @description
#' Print method for displaying summary information from distance-based bandwidth selection in 2D local polynomial regression discontinuity (RD) designs, as produced by \code{\link{rdbw2d.dist}}.
#'
#' @param x An object of class \code{rdbw2d.dist}, returned by \code{\link{rdbw2d.dist}}.
#' @param ... Additional arguments passed to the method (currently ignored).
#'
#' @return No return value. This function is called for its side effects: it prints summary information of \code{\link{rdbw2d.dist}}.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso \code{\link{rdbw2d.dist}} for distance-based bandwidth selection in 2D local polynomial RD design.
#'
#' Supported methods: \code{\link{print.rdbw2d.dist}}, \code{\link{summary.rdbw2d.dist}}.
#' @export

print.rdbw2d.dist <- function(x,...){
  cat("Call: rdbw2d.dist\n\n")

  cat(sprintf("Number of Obs.         %-10s\n", x$opt$N))
  cat(sprintf("BW type                %s\n", paste(x$opt$bwselect, "rot", sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), "rad", sep = "-")))
  cat(sprintf("Kink                   %-10s\n", x$opt$kink))
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  cat("\n")

  cat(sprintf("Number of Obs.         %-10s%-10s\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Estimand (deriv)       %-10d%-10d\n", 0, 0))
  cat(sprintf("Order est. (p)         %-10s%-10s\n", x$opt$p, x$opt$p))
  cat("\n")
  #cat("Use summary(...) to show bandwidths.\n")
}

################################################################################
#' Summary Method for Bandwidth Selection in 2D Local Polynomial RD Design (Distance-Based)
#'
#' @description
#' Summarizes bandwidth selection results from a 2D local polynomial regression discontinuity (RD) design using distance-based methods, as returned by \code{\link{rdbw2d.dist}}.
#'
#' @param object An object of class \code{rdbw2d.dist}, returned by \code{\link{rdbw2d.dist}}.
#' @param ... Optional arguments. Supported options include:
#'   \itemize{
#'     \item \code{subset}: Integer vector of indices of evaluation points to display. Defaults to all evaluation points.
#'     \item \code{sep}: Integer vector of length two. Controls spacing in the output.
#'       \code{sep[1]} controls spacing for the columns of evaluation points in the table.
#'       \code{sep[2]} controls spacing for the columns of bandwidths in the table.
#'       Default is \code{c(8, 14)}.
#'   }
#'
#' @return No return value. This function is called for its side effects: it prints a formatted summary of \code{\link{rdbw2d.dist}} results.
#'
#' @author
#' Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu} \cr
#' Rocío Titiunik, Princeton University. \email{titiunik@princeton.edu} \cr
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}
#'
#' @seealso \code{\link{rdbw2d.dist}} for bandwidth selection using 2D local polynomial RD design with distance-based methods.
#'
#' Supported methods: \code{\link{print.rdbw2d.dist}}, \code{\link{summary.rdbw2d.dist}}.
#'
#' @export

summary.rdbw2d.dist <- function(object,...) {

  x <- object
  
  args <- list(...)
  
  if (is.null(args[['subset']])) {
    subset <- NULL
  } else {
    subset <- args[['subset']]
  }

  if (is.null(args[['sep']])) {
    sep <- c(8, 14)
  } else {
    sep <- args[['sep']]
  }

  cat("Call: rdbw2d.dist\n\n")

  cat(sprintf("Number of Obs.         %-10s\n", x$opt$N))
  cat(sprintf("BW type                %s\n", paste(x$opt$bwselect, "rot", sep = "-")))
  cat(sprintf("Kernel                 %s\n", paste(tolower(x$opt$kernel), "rad", sep = "-")))
  cat(sprintf("Kink                   %-10s\n", x$opt$kink))
  cat(sprintf("VCE method             %s\n", paste(x$opt$vce, ifelse(x$opt$clustered, "-clustered", ""),sep = "")))
  cat(sprintf("Masspoints             %s\n", x$opt$masspoints))
  cat("\n")

  cat(sprintf("Number of Obs.         %-10d   %-10d\n", x$opt$N.0, x$opt$N.1))
  cat(sprintf("Estimand (deriv)       %-10d   %-10d\n", 0, 0))
  cat(sprintf("Order est. (p)         %-10d   %-10d\n", x$opt$p, x$opt$p))
  cat("\n")

  cat("Bandwidth Selection","\n")

  # Define column headers and their widths

  eval.specified <- !all(is.na(x$opt$b)) # TRUE is any entry is specified

  if (eval.specified){
    headers <- c("ID", "b1", "b2", "h0", "h1")
    col_widths <- c(4, sep[1], sep[1], sep[2], sep[2])
  } else {
    headers <- c("ID", "h0", "h1")
    col_widths <- c(4, sep[2], sep[2])
  }

  col_widths_default <- c(4, sep[1], sep[1], sep[2], sep[2])

  cat(strrep("=", sum(col_widths)), "\n")

  # Format headers using formatC for right alignment
  if (eval.specified){
    group_headers <- c(
      formatC("        Bdy Points", width = col_widths[1] + col_widths[2] + col_widths[3], format = "s", flag = "-"),
      formatC("BW Control", width = col_widths_default[4], format = "s"),
      formatC("BW Treatment", width = col_widths_default[5], format = "s")
    )
  } else{
    group_headers <- c(
      formatC("", width = col_widths[1], format = "s"),
      formatC("BW Control", width = col_widths_default[4], format = "s"),
      formatC("BW Treatment", width = col_widths_default[5], format = "s")
    )
  }

  cat(paste(group_headers, collapse = ""), "\n")

  formatted_headers <- mapply(function(h, w) formatC(h, width = w, format = "s"), headers, col_widths)
  cat(paste(formatted_headers, collapse = ""), "\n")

  # Print separator line
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
    index <- formatC(j, width = col_widths_default[1], format = "d")
    bdy1 <- ifelse(is.na(x$bws[j, "b1"]),
                   formatC("NA", width = col_widths_default[2], format = "s"),
                   formatC(x$bws[j, "b1"], format = "f", digits = 3, width = col_widths_default[2]))
    bdy2 <- ifelse(is.na(x$bws[j, "b2"]),
                   formatC("NA", width = col_widths_default[3], format = "s"),
                   formatC(x$bws[j, "b2"], format = "f", digits = 3, width = col_widths_default[3]))
    control <- ifelse(is.na(x$bws[j, "h0"]),
                      formatC("NA", width = col_widths_default[4], format = "s"),
                      formatC(x$bws[j, "h0"], format = "f", digits = 3, width = col_widths_default[4]))
    treatment <- ifelse(is.na(x$bws[j, "h1"]),
                        formatC("NA", width = col_widths_default[5], format = "s"),
                        formatC(x$bws[j, "h1"], format = "f", digits = 3, width = col_widths_default[5]))

    # Combine formatted values and print the row
    if (eval.specified){
      row_vals <- c(index, bdy1, bdy2, control, treatment)
    } else {
      row_vals <- c(index, control, treatment)
    }

    if (j %in% subset) cat(paste(row_vals, collapse = ""), "\n")
  }

  # Print closing separator line
  cat(strrep("=", sum(col_widths)), "\n")
}
