#' @title Plot scaffold
#'
#' @description Generate plots for a specific scaffold and/or region from the results of PSASS.
#'
#' @param contig_lengths_file_path Path to a contig lengths file.
#'
#' @param chromosomes_names_file_path Path to a contig names file (default NULL).
#'
#' @param prefix Prefix (including full path) to a complete dataset. If prefix is specified, it will be override individual file specifications
#' such as "window_fst_file_path" (default NULL).
#'
#' @param window_fst_file_path Path to a FST window output file (default NULL).
#'
#' @param position_fst_file_path Path to a FST positions output file (default NULL).
#'
#' @param window_snps_file_path Path to a SNPs window output file (default NULL).
#'
#' @param position_snps_file_path Path to a SNPs positions output file (default NULL).
#'
#' @param depth_file_path Path to a depth output file (default NULL).
#'
#' @param scaffold Name of the scaffold to be plotted. The name can be from the dataset (e.g. "NC_0245486.2") or from the
#' 'chromosomes_names' file (e.g. "LG4").
#'
#' @param region A vector specifying the boundaries of the region to be plotted, e.g. c(125000, 250000). If NULL, the entire scaffold
#' will be plotted (default: NULL).
#'
#' @param output.file Path to an output file in PNG format. If NULL, the plot will be drawn in the default graphic device (default: NULL).
#'
#' @param width Width of the output file if specified, in inches (default: 12).
#'
#' @param height Height of the output file if specified, in inches (default: 4 * number of plots).
#'
#' @param dpi Resolution of the output file if specified, in dpi (default: 300).
#'
#' @param tracks Tracks to be plotted. Possible values are "position_fst", "window_fst", "window_snp_males",
#' "window_snp_females", "window_snp_combined", "window_snp_ratio", "depth_males", "depth_females", "depth_combined", "depth_ratio"
#' (default: c("window_fst", "window_snp_combined", "depth_ratio")).
#'
#' @param point.size Size of the points in the plot (default: 0.1).
#'
#' @param depth.type Type of depth to be plotted, either "absolute" or "relative" (default: "absolute").
#'
#' @param depth.ratio.lines If TRUE, lines will be drawn for ratios of 2, 3/2, 2/3, and 1/2 (default: FALSE).
#'
#' @param min.depth Minimum depth to compute depth ratio.
#' The log of ratio for positions with depth lower than this value in either sex will be 0 (default: 10).
#'
#' @param major.lines.y If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param major.lines.x If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param alpha.value Transparency values for combined tracks (default: 0.25).
#'
#' @param fst.window.color Color of the FST window track (default: "grey50").
#'
#' @param males.color Color of the male window tracks (default: "dodgerblue3").
#'
#' @param females.color Color of the female window tracks (default: "firebrick2").
#'
#' @examples
#'


plot_scaffold <- function(scaffold, region = NULL,
                          contig_lengths_file_path = NULL, prefix = NULL,
                          window_fst_file_path = NULL, position_fst_file_path = NULL,
                          window_snps_file_path = NULL, position_snps_file_path = NULL,
                          depth_file_path = NULL, chromosomes_names_file_path = NULL,
                          output.file = NULL, width = 12, height = 4, dpi = 300,
                          tracks = c("window_fst", "window_snp_combined", "depth_ratio"),
                          point.size = 0.5, depth.type = "relative", min.depth = 10, depth.ratio.lines = FALSE,
                          major.lines.y = TRUE, major.lines.x = FALSE,
                          fst.window.color = "grey50", males.color = "dodgerblue3", females.color = "firebrick2", alpha.value = 0.25) {


    # Load the data
    data <- load_data_files(contig_lengths_file_path,
                            prefix = prefix,
                            window_fst_file_path = window_fst_file_path,
                            position_fst_file_path = position_fst_file_path,
                            window_snps_file_path = window_snps_file_path,
                            position_snps_file_path = position_snps_file_path,
                            depth_file_path = depth_file_path,
                            chromosomes_names_file_path = chromosomes_names_file_path,
                            plot.unplaced = TRUE)

    # Draw the plot
    draw_scaffold_plot(data,
                       scaffold = scaffold,
                       region = region,
                       output.file = output.file,
                       width = width,
                       height = height,
                       dpi = dpi,
                       tracks = tracks,
                       point.size = point.size,
                       depth.type = depth.type,
                       min.depth = min.depth,
                       depth.ratio.lines = depth.ratio.lines,
                       major.lines.y = major.lines.y,
                       major.lines.x = major.lines.x,
                       fst.window.color = fst.window.color,
                       males.color = males.color,
                       females.color = females.color,
                       alpha.value = alpha.value)
}
