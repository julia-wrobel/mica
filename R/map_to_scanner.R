#' Map niftis from one scanner to another
#'
#' This function reads in filepaths for a set of niftis from two scanners and maps niftis from one scanner to the
#' other.
#'
#' @param inpaths List of paths to where nifti objects are stored.
#' @param outpath Directory where mica normalized nifti objects will be stored.
#' @param subjects List that denotes the subject id for each image. Alignment will occur at the subject level.
#' @param scan_ids List of unique identifiers for each image.
#' @param scanner Scanner on which each image was collected. Should be character vector same length as suubjects.
#' @param map_from Character, name of scanner from which images will be mapped.
#' @param map_to Character, name of scanner to which images will be mapped.
#' @param grid_length Length of downsampled CDFs to be aligned.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @import dplyr
#' @importFrom purrr map map2 pmap
#' @importFrom tidyr nest unnest gather pivot_wider
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#' @export

map_to_scanner <- function(inpaths, outpath, subjects, scan_ids, scanner,
                           map_from, map_to,
                           grid_length = 1000, ...){

  # can only warp from one scanner to another - not from multiple scanners to another
  # test that map to and map from variables are same as in scanner variable

  subj_scan_scanner = paste(subjects, scan_ids, scanner, sep = "_")

  intensity_df = make_intensity_df(inpaths, subj_scan_scanner)

  intensity_extrema = intensity_df %>%
    filter(scanner == map_to) %>%
    group_by(subject) %>%
    summarize(intensity_minimum = min(intensity),
              intensity_maximum = max(intensity)) %>%
    ungroup()

  intensity_df = left_join(intensity_df, intensity_extrema)

  cdfs = estimate_cdf(intensity_df,
                      grid_length = grid_length)

  intensity_df = cdfs$intensity_df

  intensity_df_short = tibble(subject = rep(cdfs$intensity_df$subject, each = grid_length),
                              scan_id = rep(cdfs$intensity_df$scan_id, each = grid_length),
                              scanner = rep(cdfs$intensity_df$scanner, each = grid_length),
                              cdf = as.vector(cdfs$cdf_mat),
                              intensity = as.vector(cdfs$intensity_mat))

  intensity_df_short = intensity_df_short %>%
    nest(data = c(cdf, intensity)) %>%
    group_by(subject) %>%
    mutate(n_warps = n_distinct(scan_id) - 1) %>%
    ungroup()


  # estimate warping
  total_warps = intensity_df_short %>%
    select(subject, n_warps) %>%
    distinct() %>%
    summarize(total_warps = sum(n_warps)) %>% pull(total_warps)

  num_subjects = length(unique(subjects))

  h_inv = matrix(cdfs$intensity_mat[,1], nrow = grid_length, ncol = 1)
  colnames(h_inv) = "intensity"
  uniq_subjects = unique(subjects)

  for(i in 1:num_subjects){
    target = intensity_df_short %>% filter(subject == uniq_subjects[i], scanner == map_to)
    sources = intensity_df_short %>% filter(subject == uniq_subjects[i], scanner != map_to)

    for(scan in 1:target$n_warps){
      tstar = sources$data[[scan]]$intensity ## need to rescale intensity here?
      q_from = quantile(sources$data[[scan]]$cdf, probs = seq(0, 1, length.out = grid_length))

      scan_id = sources[scan,]$scan_id
      h_inv = cbind(h_inv, approx(target$data[[1]]$cdf, tstar, xout = q_from, rule = 2)$y)
      colnames(h_inv)[dim(h_inv)[2]] = scan_id
    }
  }

  h_inv = h_inv %>% as_tibble() %>%
    gather(scan_id, h_inv, !intensity) %>%
    select(scan_id, h_inv)

  hinv_df = intensity_df_short %>%
    filter(scanner == map_from) %>%
    unnest(cols = c(data)) %>%
    mutate(h_inv = h_inv$h_inv) %>%
    nest(data = c(cdf, intensity, h_inv)) %>% select(-n_warps)

  # upsample inverse warping functions
  intensity_df = rbind(intensity_df %>%
    filter(scanner != map_to) %>%
    arrange(scan_id, scanner) %>%
    mutate(short_data = hinv_df$data,
           data = map2(data, short_data, upsample_hinv)) %>%
    select(-short_data),
    (intensity_df %>%
    filter(scanner == map_to) %>%
    unnest(data) %>%
    mutate(h_inv = NA) %>%
    nest(data = c(intensity, voxel_position, cdf, h_inv)) )
    ) %>%
    arrange(scan_id, scanner)

  # last step is normalizing the niftis themselves
  map_from_inpaths = inpaths[which(scanner != map_to)]
  map_from_subjects = subjects[which(scanner != map_to)]
  map_from_scanids = scan_ids[which(scanner != map_to)]
  map_from_df = intensity_df %>% filter(scanner != map_to)

  hinv_ls = list(as.list(map_from_inpaths), as.list(map_from_subjects), as.list(map_from_scanids), map_from_df$data)
  image_norm = pmap(hinv_ls, .f = normalize_image, outpath = outpath)

  return(list(long_data = intensity_df, short_data = intensity_df_short, short_hinv = hinv_df))

} # end function
