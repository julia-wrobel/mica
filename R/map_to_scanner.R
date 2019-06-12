#' Map NIfTIs from one scanner to another
#'
#' This function reads in filepaths for a set of NIfTIs from two scanners and maps NIfTIs from one scanner to the
#' other.
#'
#' @param inpaths List of paths to where NIfTI objects are stored.
#' @param outpath Character, directory where mica normalized NIfTI objects will be stored.
#' @param ids List of unique identifiers for each subject. Should be character vector same length as inpaths.
#' @param scanner Scanner on which each image was collected. Should be character vector same length as inpaths.
#' @param map_from Character, name of scanner from which images will be mapped.
#' @param map_to Character, name of scanner to which images will be mapped.
#' @param white_stripe If \code{TRUE} image will be white stripe normalized
#' before it is vectorized.
#' @param type If white_stripe = TRUE, user must specify the type of image, from options
#' \code{type = c("T1", "T2", "FA", "MD", "first", "last", "largest")}.
#' @param grid_length Numeric, length of downsampled CDFs to be aligned via \code{fdasrvf::time_warping()}
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

map_to_scanner <- function(inpaths, outpath, ids, scanner, map_from, map_to,
                           white_stripe = FALSE, type = NULL, grid_length = 100, ...){

  # can only warp from one scanner to another - not from multiple scanners to another
  # test that map to and map from variables are same as in scanner variable

# check that inpaths, ids, and scanner are same length
  if(length(inpaths)!=length(ids) | length(inpaths)!=length(scanner)){
    stop("inpaths, ids, and scanner must be character vectors of same length")
  }

  num_subjs = length(unique(ids))

  if(white_stripe){
    if(is.null(type)){
      stop("Input type of image to whitestripe")
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

  intensity_maxima = intensity_df %>%
    group_by(id, site) %>%
    summarize(intensity_max = max(intensity)) %>%
    filter(site == map_to) %>% select(-site) %>%
    pull(intensity_max)
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

  intensity_df_short = intensity_df_short %>% nest(-id, -site) %>%
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

  # identify mapping labels correctly
  which_map_to = which(sort(unique(scanner)) == map_to)

  intensity_df_short = intensity_df_short %>% unnest()
  intensity_df_short[[map_to]] = if(which_map_to == 2) {
    intensity_df_short$cdf1}else{intensity_df_short$cdf}
  intensity_df_short[[map_from]] = if(which_map_to == 2) {
    intensity_df_short$cdf}else{intensity_df_short$cdf1}

  # calculate inverse warping functions
  intensity_df_short = intensity_df_short %>%
    select(-cdf, -cdf1, -intensity1) %>% # intensities are the same for both scanners
    mutate(gam = as.vector(gam)) %>%
    nest(-id) %>%
    mutate(data = map(data, inverse_warps)) %>%
    unnest() %>%
    gather(site, cdf, map_to:map_from) %>%
    nest(-id, -site) %>% arrange(id, site)

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
      nest(-id, -site, -scan)

    intensity_df = intensity_df %>% nest(-id, -site, -scan)
  }

  # last step is normalizing the NIfTIs themselves
  map_from_inpaths = inpaths[which(scanner != map_to)]
  map_from_ids = unique(ids)
  map_from_df = intensity_df %>% filter(site != map_to)

  hinv_ls = list(as.list(map_from_inpaths), as.list(map_from_ids), map_from_df$data)
  image_norm = pmap(hinv_ls, .f = normalize_image, outpath = outpath)

  return(list(long_data = intensity_df, short_data = intensity_df_short))

} # end function
