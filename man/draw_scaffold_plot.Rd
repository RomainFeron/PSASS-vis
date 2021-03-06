% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_scaffold_plot.R
\name{draw_scaffold_plot}
\alias{draw_scaffold_plot}
\title{Draw scaffold plot}
\usage{
draw_scaffold_plot(
  data,
  scaffold,
  region = NULL,
  output.file = NULL,
  width = 12,
  height = 4,
  dpi = 300,
  tracks = c("window_fst", "window_snp_males", "window_snp_females", "depth_ratio"),
  point.size = 0.5,
  depth.type = "absolute",
  min.depth = 10,
  depth.ratio.lines = FALSE,
  major.lines.y = TRUE,
  major.lines.x = FALSE,
  fst.window.color = "grey50",
  males.color = "dodgerblue3",
  females.color = "firebrick2",
  alpha.value = 0.25
)
}
\arguments{
\item{data}{A PSASS data structure obtained with the \code{\link{load_data_files}} function.}

\item{scaffold}{Name of the scaffold to be plotted. The name can be from the dataset (e.g. "NC_0245486.2") or from the
'chromosomes_names' file (e.g. "LG4").}

\item{region}{A vector specifying the boundaries of the region to be plotted, e.g. c(125000, 250000). If NULL, the entire scaffold
will be plotted (default: NULL).}

\item{output.file}{Path to an output file in PNG format. If NULL, the plot will be drawn in the default graphic device (default: NULL).}

\item{width}{Width of the output file if specified, in inches (default: 12).}

\item{height}{Height of the output file if specified, in inches (default: 4 * number of plots).}

\item{dpi}{Resolution of the output file if specified, in dpi (default: 300).}

\item{tracks}{Tracks to be plotted. Possible values are "position_fst", "window_fst", "window_snp_males",
"window_snp_females", "window_snp_combined", "window_snp_ratio", "depth_males", "depth_females", "depth_combined", "depth_ratio"
(default: c("window_fst", "window_snp_males", "window_snp_females", "depth_ratio")).}

\item{point.size}{Size of the points in the plot (default: 0.1).}

\item{depth.type}{Type of depth to be plotted, either "absolute" or "relative" (default: "absolute").}

\item{min.depth}{Minimum depth to compute depth ratio.
The log of ratio for positions with depth lower than this value in either sex will be 0 (default: 10).}

\item{depth.ratio.lines}{If TRUE, lines will be drawn for ratios of 2, 3/2, 2/3, and 1/2 (default: FALSE).}

\item{major.lines.y}{If TRUE, major grid lines will be plotted for the y axis (default: TRUE).}

\item{major.lines.x}{If TRUE, major grid lines will be plotted for the y axis (default: TRUE).}

\item{fst.window.color}{Color of the FST window track (default: "grey50").}

\item{males.color}{Color of the male window tracks (default: "dodgerblue3").}

\item{females.color}{Color of the female window tracks (default: "firebrick2").}

\item{alpha.value}{Transparency values for combined tracks (default: 0.25).}
}
\description{
Draw a scaffold plot from the PSASS data
}
\examples{

}
