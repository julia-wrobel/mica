#' Calulate eCDF and downsample to regular grid.
#'
#' Function used to calculate empirical CDF of images and downsample to a
#' grid of regularly-spaced intensity values.
#'
#' @param intensity_df Dataframe of vectorized intensity values.
#' @param grid_length Length of downsampled CDFs to be aligned
#' @param ... Additional arguments passed to or from other functions.
#'
#' @import dplyr
#' @importFrom tidyr nest unnest
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#'
#' @export
estimate_cdf <- function(intensity_df,
                         grid_length = 1000, ...){

  intensity_df = intensity_df %>%
    nest(data = c(intensity, voxel_position)) %>%
    mutate(data = map(data, calculate_cdf))

  # downsample cdf to smaller, regular grid
  cdf_mat = intensity_mat = matrix(0, nrow = grid_length, ncol = dim(intensity_df)[1])

  for(i in 1:dim(intensity_df)[1]){

    intensity_minimum = intensity_df$intensity_minimum[i]
    intensity_maximum = intensity_df$intensity_maximum[i]
    intensity_grid = seq(intensity_minimum, intensity_maximum, length.out = grid_length)

    cdf_mat[, i] = approx(intensity_df$data[[i]]$intensity,
                          intensity_df$data[[i]]$cdf,
                          xout = intensity_grid, rule = 2, ties = mean)$y

    intensity_mat[, i] = intensity_grid
  }
  ls = list(intensity_df = intensity_df, cdf_mat = cdf_mat, intensity_mat = intensity_mat)

  return(ls)
}
