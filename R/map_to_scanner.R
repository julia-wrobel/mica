#' Map niftis from one scanner to another
#'
#' This function reads in filepaths for a set of niftis from two scanners and maps niftis from one scanner to the
#' other.
#'
#' @param inpaths List of paths to where nifti objects are stored.
#' @param outpath Directory where mica normalized nifti objects will be stored.
#' @param ids List of unique identifiers for each subject.
#' @param scanner Scanner on which each image was collected. Should be character vector same length as ids.
#' @param map_from Character, name of scanner from which images will be mapped.
#' @param map_to Character, name of scanner to which images will be mapped.
#' @param white_stripe If \code{TRUE} image will be white stripe normalized
#' before it is vectorized.
#' @param type If white_stripe = TRUE, user must specify the type of image, from options
#' \code{type = c("T1", "T2", "FA", "MD", "first", "last", "largest")}.
#' @param grid_length Length of downsampled CDFs to be aligned via \code{fdasrvf::time_warping()}
#' @param ... Additional arguments passed to or from other functions.
#'
#' @import dplyr
#' @importFrom purrr map
#' @importFrom fdasrvf pair_align_functions
#' @importFrom tidyr nest unnest gather spread
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export

map_to_scanner <- function(inpaths, outpath, ids, scanner, map_from, map_to,
                           white_stripe = FALSE, type = NULL, grid_length = 100, ...){

  # can only warp from one scanner to another - not from multiple scanners to another
  # test that map to and map from variables are same as in scanner variable

  num_subjs = length(unique(ids))

  if(white_stripe){
    if(is.null(type)){
      stop("Input type of image to whitestripe.")
    }
  }

  intensity_df = make_intensity_df(inpaths, ids, white_stripe = white_stripe,
                                   type = type, sitenames = scanner)

  intensity_maxima = intensity_df %>%
    filter(site == map_to) %>% group_by(id) %>%
    summarize(intensity_max = max(intensity)) %>% ungroup() %>% pull(intensity_max)

  intensity_maxima = rep(intensity_maxima, each = 2)

  cdfs = estimate_cdf(intensity_df, intensity_maximum = intensity_maxima,
                      rescale_intensities = FALSE,
                      white_stripe = white_stripe, type = type,
                      grid_length = grid_length)

  intensity_df = cdfs$intensity_df

  intensity_df_short = tibble(id = rep(cdfs$intensity_df$id, each = grid_length),
                            site = rep(cdfs$intensity_df$site, each = grid_length),
                            cdf = as.vector(cdfs$cdf_mat),
                            intensity = as.vector(cdfs$intensity_mat))

  intensity_df_short = intensity_df_short %>% nest(cdf, intensity) %>%
    spread(site, data)

  # estimate warping
  gam = matrix(NA, nrow = grid_length, ncol = num_subjs)

  # warping old scan to new scan
  for(i in 1:num_subjs){
    id = unique(intensity_df$id)[i]
    warp_obj = pair_align_functions(intensity_df_short[[map_to]][[i]]$cdf,
                                    intensity_df_short[[map_from]][[i]]$cdf,
                                    intensity_df_short[[map_to]][[i]]$intensity)
    gam[,i] = warp_obj$gam
  }

  which_map_to = which(sort(unique(scanner)) == map_to)

  intensity_df_short = intensity_df_short %>% unnest()
  intensity_df_short[[map_to]] = if(which_map_to == 2) {
    intensity_df_short$cdf1}else{intensity_df_short$cdf}
  intensity_df_short[[map_from]] = if(which_map_to == 2) {
    intensity_df_short$cdf}else{intensity_df_short$cdf1}

  intensity_df_short = intensity_df_short %>%
    select(-cdf, -cdf1, -intensity1) %>% # intensities are the same for both scanners
    mutate(gam = as.vector(gam)) %>%
    nest(-id) %>%
    mutate(data = map(data, inverse_warps)) %>%
    unnest() %>%
    gather(scanner, cdf, map_to:map_from)

  intensity_df_short

} # end function