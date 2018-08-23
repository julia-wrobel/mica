#' create labeled dataframe of vectorized voxel intensities for a list of nifti objects
#'
#' Function used along with \code{mica::vectorize_image()} to vectorize images and place in a
#' dataframe, where names for each image are drawn from the filepath of the
#' image.
#'
#' @param inpaths List of paths to where nifti objects are stored.
#' @param outpath Directory where mica normalized nifti objects will be stored
#' @param ids List of unique identifiers for each image.
#' @param intensity_maximum Maximum value of intensity for creating grid over which to evaluate CDF and
#' estimate warping functions. Needs to be chosen based on inspection of data.
#' @param rescale_intensities If \code{TRUE}, intensities will be rescaled by their 99.9 percent quantile.
#' @param white_stripe If \code{TRUE} image will be white stripe normalized
#' before it is vectorized.
#' @param type If white_stripe = TRUE, user must specify the type of image, from options
#' \code{type = c("T1", "T2", "FA", "MD", "first", "last", "largest")}.
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
##' ids = paste(sitenames, scan_nos, sep = "_")
##'
##' intensity_data = make_intensities_df(filenames, ids, sitenames, scan_nos)
##' }
#'
#' @return a data frame of vectorized images with columns \code{intensity}, \code{id},
#' \code{site}, and \code{scan}, \code{cdf}.
#' @export

map_to_mean <- function(inpaths, outpath, ids, intensity_maximum, rescale_intensities = FALSE,
                        white_stripe = FALSE, type = NULL, grid_length = 100, ...){

  if(white_stripe){
    if(is.null(type)){
      stop("Input type of image to whitestripe.")
    }
  }

  intensity_df = make_intensity_df(inpaths, ids, white_stripe = white_stripe, type = type)

  intensity_grid = seq(0, intensity_maximum, length.out = grid_length)
  cdfs = estimate_cdf(intensity_df, intensity_maximum, rescale_intensities = rescale_intensities,
                      white_stripe = white_stripe, type = type,
                      grid_length = grid_length)

  intensity_df = cdfs$intensity_df

  # estimate warping
  srvf_obj = suppressMessages(time_warping(cdfs$cdf_mat, time = intensity_grid, showplot = FALSE))

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

  if(white_stripe){
    intensity_df = intensity_df %>%
      unnest(data) %>%
      mutate(h_inv = h_inv + min(intensity_ws))

    intensity_df_short = intensity_df_short %>%
      unnest(data) %>%
      mutate(intensity_ws = intensity + min(intensity_df$intensity_ws),
             h_inv = h_inv + min(intensity_df$intensity_ws)) %>%
      nest(-id, -site, -scan)

    intensity_df = intensity_df %>% nest(-id, -site, -scan)

  }

  # last step is normalizing the niftis themselves
  hinv_ls = list(as.list(inpaths), as.list(ids), intensity_df$data)
  image_norm = pmap(hinv_ls, .f = normalize_image, outpath = outpath)

  return(list(long_data = intensity_df, short_data = intensity_df_short))
}
