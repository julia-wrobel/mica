#' Normalize a set of niftis to their Karcher mean.
#'
#' This function reads in filepaths for a set of niftis and outputs mica-normalized niftis to a specified folder.
#' Also returns dataframes of vectorized niftis and warping functions.
#'
#' @param inpaths List of paths to where nifti objects are stored.
#' @param outpath Directory where mica normalized nifti objects will be stored
#' @param subjects List that denotes the subject id for each image. Alignment will occur at the subject level.
#' @param scan_id List of unique identifiers for each image.
#' @param intensity_maximum Maximum value of intensity for creating grid over which to evaluate CDF and
#' estimate warping functions. Needs to be chosen based on inspection of data.
#' @param rescale_intensities If \code{TRUE}, intensities will be rescaled by their 99.9 percent quantile.
#' @param grid_length Length of downsampled CDFs to be aligned via \code{fdasrvf::time_warping()}
#' @param ... Additional arguments passed to or from other functions.
#'
#' @import dplyr
#' @importFrom purrr map2 map_df map pmap
#' @importFrom fdasrvf time_warping
#' @importFrom tidyr nest unnest
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @examples
##' \dontrun{
##' filenames = list.files(path = "../Documents/", pattern = "flair_reg", recursive = TRUE)
##' filenames = paste0("../Documents/", filenames)
##' sitenames = lapply(filenames, function(f) substring(f, 64, 66))
##' scan_nos = lapply(filenames, function(f) substring(str_split(f, "Scan_")[[1]][2], 1, 1))
##' subjects = paste(sitenames, scan_nos, sep = "_")
##'
##' intensity_data = make_intensities_df(filenames, subjects, sitenames, scan_nos)
##' }
#'
#' @return a data frame of vectorized images with columns \code{intensity}, \code{id},
#' \code{site}, and \code{scan}, \code{cdf}.
#' @export

map_to_mean <- function(inpaths, outpath, subjects, intensity_maximum, rescale_intensities = FALSE,
                        grid_length = 100, ...){

  if(white_stripe){
    if(is.null(type)){
      stop("Input type of image to whitestripe.")
    }
  }

  intensity_df = make_intensity_df(inpaths, subjects, white_stripe = white_stripe, type = type)

  intensity_grid = seq(0, intensity_maximum, length.out = grid_length)
  cdfs = estimate_cdf(intensity_df, intensity_maximum, rescale_intensities = rescale_intensities,
                      white_stripe = white_stripe, type = type,
                      grid_length = grid_length)

  intensity_df = cdfs$intensity_df

  # estimate warping
  srvf_obj = suppressMessages(time_warping(cdfs$cdf_mat, time = intensity_grid,
                                           showplot = TRUE))

  intensity_df_short = tibble(
    id = rep(cdfs$intensity_df$id, each = grid_length),
    site = rep(cdfs$intensity_df$site, each = grid_length),
    scan = rep(cdfs$intensity_df$scan, each = grid_length),
    intensity = rep(intensity_grid, dim(cdfs$cdf_mat)[2]),
    cdf = as.vector(cdfs$cdf_mat),
    gam = as.vector(srvf_obj$gam)
  )

  # calculate inverse warping functions
  intensity_df_short = intensity_df_short %>%
    nest(intensity, cdf, gam) %>%
    mutate(data = map(data, inverse_warps))

  # upsample inverse warping functions
  intensity_df = intensity_df %>%
    mutate(short_data = intensity_df_short$data,
           data = map2(data, short_data, upsample_hinv)) %>%
    select(-short_data)


  # last step is normalizing the niftis themselves
  hinv_ls = list(as.list(inpaths), as.list(subjects), intensity_df$data)
  image_norm = pmap(hinv_ls, .f = normalize_image, outpath = outpath)

  return(list(long_data = intensity_df, short_data = intensity_df_short))
}
