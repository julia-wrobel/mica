#' create labeled dataframe of vectorized voxel intensities for a list of nifti objects
#'
#' Function used along with \code{mica::vectorize_image()} to vectorize images and place in a
#' dataframe, where names for each image are drawn from the filepath of the
#' image.
#'
#' @param filepaths List of paths to where nifti objects are stored.
#' @param subj_scan_scanner Character valued ID for subject, scan number, and scanner, in the format
#' subj_scan_scanner.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom purrr map2 map_df
#' @importFrom tibble as_tibble
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
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

make_intensity_df <- function(filepaths, subj_scan_scanner,...){

  intensities = map2_dfr(filepaths, subj_scan_scanner, vectorize_image)

  as_tibble(intensities)
}
