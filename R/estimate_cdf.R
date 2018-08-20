

## next we want to generate the actual cdfs based on a data frame that is fed in with the intensities.
# How do I want to select where to warp the intensities to? If we're whitestriping first than may not need
# the scaling step. Maybe just warp to range where most intensities are at their 99th percentile?
# look at some pictures to determine this. For naims data we already have a set range,
# though after whitestriping that range will be different.

# then:
  # calculate CDF (should I do this is prior function? do we need full length cdf?)


#' create labeled dataframe of vectorized voxel intensities for a list of nifti objects
#'
#' Function used along with \code{mica::vectorize_image()} to vectorize images and place in a
#' dataframe, where names for each image are drawn from the filepath of the
#' image.
#'
#' @param intensity_df Dataframe of vectorized intensity values.
#' @param intensity_maximum Maximum value of intensity for creating grid over which to evaluate CDF and
#' estimate warping functions. Needs to be chosen based on inspection of data.
#' @param rescale_intensities If \code{TRUE}, intensities will be rescaled by their 99.9% quantile.
#' Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @import dplyr
#' @importFrom tidyr nest unnest
#' @importFrom purrr map
#' @importFrom magrittr %>%
#'
estimate_cdf <- function(intensity_df, intensity_maximum, rescale_intensities = FALSE, ...){

  intensity_grid = seq(0, intensity_maximum, length.out = 100)

  if(rescale_intensities){ # if not whitestriping or if intensities are very different across images, rescale intensities so warping functions are identifiable
    intensity_df2 = intensity_df %>%
      group_by(id) %>%
      mutate(quantile99 = quantile(intensity, .999),
             intensity_stretch = intensity * intensity_maximum / quantile99) %>%
      ungroup() %>%
      nest(intensity) %>%
      mutate(cdf = map(data, interpolate_cdf, intensity_grid = intensity_grid,
                       intensity_values = intensity_stretch, long_to_long = TRUE)) %>%
      unnest(data, cdf)
  }else if(white_stripe){
    intensity_df = intensity_df %>%
      mutate(intensity_positive = intensity + min(intensity)) %>%
      nest(intensity) %>%
      mutate(cdf = map(data, interpolate_cdf, intensity_grid = intensity_grid,
                       intensity_values = intensity_positive, long_to_long = TRUE)) %>%
      unnest(data, cdf)
  }else{
    intensity_df = intensity_df %>%
      nest(intensity) %>%
      mutate(cdf = map(data, interpolate_cdf, intensity_grid = intensity_grid,
                       intensity_values = intensity, long_to_long = TRUE)) %>%
      unnest(data, cdf)
  }

  # interpolate to smaller, even grid (did I need to cal)

} # end of overall function
