################################################################################
#' @title rd2d: Two Dimensional Local Polynomial Regression Discontinuity Design
#'
#' @description This package implements estimation and inference procedures for boundary regression discontinuity (RD) designs
#' using local polynomial methods, based on either bivariate coordinates or distance-based approaches.
#' Methods are developed in Cattaneo, Titiunik, and Yu (2025) 
#' <https://mdcattaneo.github.io/papers/Cattaneo-Titiunik-Yu_2025_BoundaryRD.pdf>.
#' 
#' Included functions are: \link{rd2d} for inference and estimation based on bivariate coordinates, 
#' \link{rdbw2d} for data-driven bandwidth selection based on bivariate coordinates,  
#' \link{rd2d.dist} for distance-based inference and estimation,
#' \link{rdbw2d.dist} for distance-based bandwidth selection.
#' 
#' \code{print()} and \code{summary()} methods are available all four functions.
#' 
#' Related Stata, R, and Python packages useful for inference in RD designs are described in the following website:
#' 
#' \href{ https://rdpackages.github.io/}{ https://rdpackages.github.io/}
#' 
#' For an introduction to regression discontinuity design, see \href{https://www.cambridge.org/core/elements/practical-introduction-to-regression-discontinuity-designs/C6A70A32359115510AAC370A7869AE2F}{Cattaneo (2024)} and references therein.
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}.
#' Ruiqi Rae Yu, Princeton University. \email{rae.yu@princeton.edu}.
#'
#' @references
#' \itemize{
#' \item{\href{https://mdcattaneo.github.io/papers/Cattaneo-Titiunik-Yu_2025_BoundaryRD.pdf}{Cattaneo, M. D., Titiunik, R., Yu, R. R. (2025).}
#' Estimation and Inference in Boundary Discontinuity Designs}
#' \item{\href{https://www.cambridge.org/core/elements/practical-introduction-to-regression-discontinuity-designs/C6A70A32359115510AAC370A7869AE2F}{Cattaneo, M. D., Idrobo, N., Titiunik, R. (2024).}
#' A Practical Introduction to Regression Discontinuity Designs: Extensions}
#' }
#'
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats D
#' @importFrom stats integrate
#' @importFrom stats optimize
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom stats sd
#' @importFrom stats as.formula
#' @importFrom stats complete.cases
#' @importFrom stats cov
#' @importFrom stats lm
#' @importFrom stats median
#' @importFrom stats predict
#' @importFrom stats density
#' @importFrom MASS mvrnorm
#' @importFrom MASS ginv
#' @importFrom expm sqrtm

#' @import ggplot2
#'
#' @aliases rd2d-package
"_PACKAGE"
