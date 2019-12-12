#' @title Plot manhattan
#'
#' @description Generate a manhattan plot for a specific track from the results of PSASS.
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
#' @param output.file Path to an output file in PNG format. If NULL, the plot will be drawn in the default graphic device (default: NULL).
#'
#' @param width Width of the output file if specified, in inches (default: 14).
#'
#' @param height Height of the output file if specified, in inches (default: 8).
#'
#' @param dpi Resolution of the output file if specified, in dpi (default: 300).
#'
#' @param track Track to be plotted. Possible values are "position_fst", "window_fst", "window_snp_males", "window_snp_females", "depth_males", "depth_females"
#' (default: "window_fst").
#'
#' @param lg.numbers If TRUE, chromosomes / LGs will be labeled with numbers instead of names to increase readability (default: FALSE).
#'
#' @param point.size Size of a point in the plot (default 0.5)
#'
#' @param point.palette Color palette for the dots (default c("dodgerblue3", "darkgoldenrod2"))
#'
#' @param background.palette Color palette for the background (default c("grey85", "grey100"))
#'
#' @param depth.type Type of depth to be plotted, either "absolute" or "relative" (default: "absolute").
#'
#' @param min.depth Minimum depth to compute depth ratio.
#' The ratio for positions with depth lower than this value in either sex will be 1 (default: 10).
#'
#' @examples
#'
#' c_length <- "genome.fasta.fai"
#' c_names <- "names.tsv"
#' prefix <- "psass"
#' plot_manhattan(c_length, prefix=prefix, chromosomes_names_file_path = c_names, lg.numbers = TRUE)

plot_manhattan <- function(contig_lengths_file_path,
                           prefix = NULL,
                           window_fst_file_path = NULL, position_fst_file_path = NULL,
                           window_snps_file_path = NULL, position_snps_file_path = NULL,
                           depth_file_path = NULL, chromosomes_names_file_path = NULL,
                           output.file = NULL, width = 14, height = 8, dpi = 300,
                           track = "window_fst", lg.numbers = FALSE,
                           point.size = 0.5, point.palette = c("dodgerblue3", "darkgoldenrod2"), background.palette = c("grey85", "grey100"),
                           depth.type = "absolute", min.depth = 10) {


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

    draw_manhattan_plot(data,
                        output.file = output.file,
                        width = width,
                        height = height,
                        dpi = dpi,
                        track = track,
                        lg.numbers = lg.numbers,
                        point.size = point.size,
                        point.palette = point.palette,
                        background.palette = background.palette,
                        depth.type = depth.type,
                        min.depth = min.depth)
}
