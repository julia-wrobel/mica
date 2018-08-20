#' Interpolate CDF
#'
#' Given a set of intensity values for an image, will interpolate the empirical CDF to a new domain.
#' Can calculate a new grid that is sparse (if long_to_short = TRUE) or denser
#' (if long_to_short = FALSE)
#'
#' @param data dataset with intensity values.
#' @param intensity_values values of intensity to determine CDF of for alignment.
#' @param intensity_grid set of intensities shorter than the length of the original data
#' @param long_to_long If TRUE, outputs CDFs for full intensity grid.
#' If FALSE, takes sparse intensity grid as input and outputs CDF the length of the original data.
#'
#' @importFrom stats ecdf
#' @export
interpolate_cdf = function(data, intensity_values, intensity_grid, long_to_long = FALSE){
  long_intensities = data[[intensity_values]]
  short_intensities = intensity_grid

  if(long_to_long){
    ecdf(long_intensities)(long_intensities)
  }else{
    ecdf(long_intensities)(short_intensities)
  }
}
