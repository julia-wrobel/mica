#' Map niftis from one scanner to another
#'
#' This function reads in filepaths for a set of niftis from two scanners and maps niftis from one scanner to the
#' other.
#'
#' @param inpaths List of paths to where nifti objects are stored.
#' @param outpath Directory where mica normalized nifti objects will be stored.
#' @param ids List of unique identifiers for each subject.
#' @param scanner Scanner on which each image was collected. Should be character vector same length as ids.
#' @param intensity_maximum Maximum value of intensity for creating grid over which to evaluate CDF and
#' estimate warping functions. Needs to be chosen based on inspection of data.
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
#' @importFrom purrr map map2 pmap
#' @importFrom fdasrvf pair_align_functions
#' @importFrom tidyr nest unnest gather spread
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' @export

map_to_scanner <- function(inpaths, outpath, ids, scanner, intensity_maximum = NULL, map_from, map_to,
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

  # add intensity_maxima to the dataset
  # intensity_maxima = intensity_df %>%
  #   group_by(id, site) %>%
  #   summarize(intensity_max = max(intensity)) %>%
  #   filter(site == map_to) %>% select(-site)

  #intensity_df  = left_join(intensity_df, intensity_maxima)
  #intensity_maxima = rep(unique(intensity_df$intensity_max), each = 2)
  if(is.null(intensity_maximum)){
    intensity_maximum = intensity_df %>%
      group_by(id, site) %>%
      summarize(intensity_max = max(intensity)) %>%
      filter(site == map_to) %>% select(-site) %>%
      pull(intensity_max)
    intensity_maximum = rep(intensity_maximum, each = 2)
  }

  cdfs = estimate_cdf(intensity_df, intensity_maximum,
                      rescale_intensities = FALSE,
                      white_stripe = white_stripe, type = type,
                      grid_length = grid_length)

  intensity_df = cdfs$intensity_df

  intensity_df_short = tibble(id = rep(cdfs$intensity_df$id, each = grid_length),
                            site = rep(cdfs$intensity_df$site, each = grid_length),
                            cdf = as.vector(cdfs$cdf_mat),
                            intensity = as.vector(cdfs$intensity_mat))

  intensity_df_short = intensity_df_short %>%
    nest(data = c(cdf, intensity)) %>%
    pivot_wider(names_from = site, values_from = data)

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

  intensity_df_short = intensity_df_short %>% unnest(cols = c(prisma)) %>%
    rename(prisma = cdf, intensity_prisma = intensity) %>%
    nest(prisma = c(prisma, intensity_prisma)) %>%
    unnest(cols = c(trio, prisma)) %>%
    rename(trio = cdf, intensity_trio = intensity)

  # calculate inverse warping functions
  intensity_df_short = intensity_df_short %>%
    select( -intensity_trio) %>% # intensities are the same for both scanners
    rename(intensity = intensity_prisma) %>%
    mutate(gam = as.vector(gam)) %>%
    nest(data = c(trio, prisma, intensity, gam)) %>%
    mutate(data = map(data, inverse_warps)) %>%
    unnest(cols = c(data)) %>%
    gather(site, cdf, map_to:map_from) %>%
    nest(data = c(intensity, gam, h_inv, cdf)) %>%
    arrange(id, site)

  # upsample inverse warping functions
  intensity_df = intensity_df %>%
    mutate(short_data = intensity_df_short$data,
           data = map2(data, short_data, upsample_hinv)) %>%
    select(-short_data) %>% arrange(id, site)

  if(white_stripe){
    intensity_df = intensity_df %>%
      unnest(data) %>%
      mutate(h_inv = h_inv + min(intensity_ws))

    intensity_df_short = intensity_df_short %>%
      unnest(data) %>%
      mutate(intensity_ws = intensity + min(intensity_df$intensity_ws),
             h_inv = h_inv + min(intensity_df$intensity_ws)) %>%
      nest(-id, -site)

    intensity_df = intensity_df %>% nest(-id, -site, -scan)
  }

  # last step is normalizing the niftis themselves
  map_from_inpaths = inpaths[which(scanner != map_to)]
  map_from_ids = unique(ids)
  map_from_df = intensity_df %>% filter(site != map_to)

  hinv_ls = list(as.list(map_from_inpaths), as.list(map_from_ids), map_from_df$data)
  image_norm = pmap(hinv_ls, .f = normalize_image, outpath = outpath)

  return(list(long_data = intensity_df, short_data = intensity_df_short))

} # end function
