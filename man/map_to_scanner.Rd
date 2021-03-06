% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_to_scanner.R
\name{map_to_scanner}
\alias{map_to_scanner}
\title{Map niftis from one scanner to another}
\usage{
map_to_scanner(
  inpaths,
  outpath,
  subjects,
  scan_ids,
  scanner,
  map_from,
  map_to,
  grid_length = 1000,
  ...
)
}
\arguments{
\item{inpaths}{List of paths to where nifti objects are stored.}

\item{outpath}{Directory where mica normalized nifti objects will be stored.}

\item{subjects}{List that denotes the subject id for each image. Alignment will occur at the subject level.}

\item{scan_ids}{List of unique identifiers for each image.}

\item{scanner}{Scanner on which each image was collected. Should be character vector same length as suubjects.}

\item{map_from}{Character, name of scanner from which images will be mapped.}

\item{map_to}{Character, name of scanner to which images will be mapped.}

\item{grid_length}{Length of downsampled CDFs to be aligned.}

\item{...}{Additional arguments passed to or from other functions.}
}
\description{
This function reads in filepaths for a set of niftis from two scanners and maps niftis from one scanner to the
other.
}
\author{
Julia Wrobel \email{julia.wrobel@cuanschutz.edu}
}
