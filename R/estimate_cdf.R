#' Calulate eCDF and downsample to regular grid.
#'
#' Function used to calculate empirical CDF of images and downsample to a
#' grid of regularly-spaced intensity values.
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
#' @importFrom stats quantile
#'
#' @export
estimate_cdf <- function(intensity_df, intensity_maximum = NULL, rescale_intensities = FALSE,
                         grid_length = 100, ...){

### automatically determine intensity_maximum

  # test that length of intensity_maximum is either 1 or length of unique ids
  if(rescale_intensities){ # if not whitestriping or if intensities are very different across images, rescale intensities so warping functions are identifiable
    intensity_df = intensity_df %>%
      group_by(id) %>%
      mutate(quantile99 = quantile(intensity, .999),
             intensity_raw = intensity,
             intensity = intensity * intensity_maximum / quantile99) %>%
      ungroup()
    # intensity_df = intensity_df %>%
    #   group_by(id) %>%
    #   mutate(intensity_raw = intensity,
    #          intensity = intensity * intensity_max / max(intensity)) %>%
    #   ungroup()
  }

  if(min(intensity_df$intensity)<0){
    intensity_df = intensity_df %>%
      mutate(intensity_ws = intensity,
             intensity = intensity - min(intensity))
  }

  intensity_df = intensity_df %>%
    nest(-id, -site, -scan) %>%
    mutate(data = map(data, calculate_cdf))

  # downsample cdf to smaller, regular grid
  cdf_mat = intensity_mat = matrix(0, nrow = grid_length, ncol = dim(intensity_df)[1])

  for(i in 1:dim(intensity_df)[1]){
    if(length(intensity_maximum) > 1){
      intensity_grid = seq(0, intensity_maximum[i], length.out = grid_length)
    }else{
      intensity_grid = seq(0, intensity_maximum, length.out = grid_length)
    }
    cdf_mat[, i] = approx(intensity_df$data[[i]]$intensity,
                          intensity_df$data[[i]]$cdf,
                          xout = intensity_grid, rule = 2)$y

    intensity_mat[, i] = intensity_grid
  }
  ls = list(intensity_df = intensity_df, cdf_mat = cdf_mat, intensity_mat = intensity_mat)

  return(ls)
}
