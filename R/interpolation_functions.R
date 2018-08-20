#' Calculate empirical CDF
#'
#' Given a set of intensity values for an image, will interpolate the
#' empirical CDF to a new domain.
#'
#' @importFrom stats ecdf
#' @export
#'
calculate_cdf <- function(data){
  intensities = data$intensity
  data$CDF = ecdf(intensities)(intensities)
  data
}
