#' Calculate empirical CDF
#'
#' Given a dataset with intensity values for an image, will interpolate the
#' empirical CDF to a new domain.
#'
#' @param data A dataframe with column \code{intensity}.
#'
#' @importFrom stats ecdf
#' @export

calculate_cdf <- function(data){
  intensities = data$intensity
  data$cdf = ecdf(intensities)(intensities)
  data
}


#' Calculate inverse warping functions
#'
#' Given a dataset with columns \code{intensity} and \code{gam}  will
#' calculate the inverse warping functions
#'
#' @param data A dataframe with columns \code{intensity} and \code{gam}.
#'
#' @importFrom stats approx
#' @export

inverse_warps = function(data){
  tstar = data[["intensity"]]
  gam = data[["gam"]]

  gamI = approx(gam*(max(tstar) - min(tstar)) + min(tstar), tstar, xout = tstar)$y

  data$h_inv = gamI
  data
}


#' Upsample inverse warping functions
#'
#' Given a dataset with columns \code{intensity} and \code{gam}  will
#' upsample the inverse warping functions
#'
#' @param long_data A dataframe with columns of nifti-length CDFs and intensity values.
#' @param short_data A downsampled dataframe with columns \code{intensity} and \code{h_inv}.
#'
#' @importFrom stats approx
#' @export
upsample_hinv = function(long_data, short_data){
  long_data$h_inv = approx(short_data$intensity, short_data$h_inv,
                           xout = long_data$intensity,
                           yleft = min(short_data$h_inv),
                           yright = max(short_data$h_inv))$y

  long_data
}
