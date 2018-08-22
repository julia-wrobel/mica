#' Calulate eCDF and downsample to regular grid
#'
#' Function used to calculate empirical CDF of images and downsample to a
#' 100-point grid of regularly-spaced intensity values.
#'
#' @param intensity_df Dataframe of vectorized intensity values.
#' @param intensity_maximum Maximum value of intensity for creating grid over which to evaluate CDF and
#' estimate warping functions. Needs to be chosen based on inspection of data.
#' @param rescale_intensities If \code{TRUE}, intensities will be rescaled by their 99.9% quantile.
#' Defaults to \code{FALSE}.
#' @param white_stripe If \code{TRUE} image will be white stripe normalized
#' @param type If \code{white_stripe = TRUE}, user must specify the type of image, from options
#' \code{type = c("T1", "T2", "FA", "MD", "first", "last", "largest")}.
#' @param grid_length Length of downsampled CDFs to be aligned via \code{fdasrvf::time_warping()}
#' @param ... Additional arguments passed to or from other functions.
#'
#' @import dplyr
#' @importFrom tidyr nest unnest
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
#' @export
estimate_cdf <- function(intensity_df, intensity_maximum, rescale_intensities = FALSE,
                         white_stripe = FALSE, type = NULL, grid_length = 100, ...){

  if(rescale_intensities){ # if not whitestriping or if intensities are very different across images, rescale intensities so warping functions are identifiable
    intensity_df = intensity_df %>%
      group_by(id) %>%
      mutate(quantile99 = quantile(intensity, .999),
             intensity_raw = intensity,
             intensity = intensity * intensity_maximum / quantile99) %>%
      ungroup()
  }

  if(white_stripe){
    if(is.null(type)){
      stop("Input type of image to whitestripe.")
    }

    intensity_df = intensity_df %>%
      mutate(intensity_ws = intensity,
             intensity = intensity - min(intensity))
  }

  intensity_df = intensity_df %>%
    nest(-id, -site, -scan) %>%
    mutate(data = map(data, calculate_cdf))

  # downsample cdf to smaller, regular grid
  cdf_mat = matrix(0, nrow = grid_length, ncol = dim(intensity_df)[1])
  intensity_grid = seq(0, intensity_maximum, length.out = grid_length)
  for(i in 1:dim(intensity_df)[1]){
    cdf_mat[, i] = approx(intensity_df$data[[i]]$intensity,
                          intensity_df$data[[i]]$CDF,
                          xout = intensity_grid, rule = 2)$y
  }
  ls = list(intensity_df = intensity_df, cdf_mat = cdf_mat)

  return(ls)
}
