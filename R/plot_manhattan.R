#' @title Manhattan plot
#'
#' @description Generate a manhattan plot with multiple tracks from a genomic data file
#'
#' @param input.file.path Path to a genomic data input file (e.g. result of PSASS or RADSex)
#'
#' @param tracks List of tracks to plot. Tracks can be generated with the \code{\link{manhattan_track}} function
#'
#' @param chromosomes.file.path Path to a tabulated file specifying chromosome names (default: NULL)
#'
#' @param detect.chromosomes If TRUE, will consider contigs starting with LG, CH, or NC as chromosomes if no chromosomes were specified (default: TRUE)
#'
#' @param unplaced.label Label for unplaced contigs (default: "Unplaced")
#'
#' @param output.file Path to an output file for the generated manhattan plot, or NULL to plot in the current R device (default: NULL)
#'
#' @param width Plot width when plotting to an output file, in inches (default: 12)
#'
#' @param track.height Height of a single track when plotting to an output file, in inches (default: 6)
#'
#' @param res Image resolution when plotting to an output file, in dpi (default: 300)
#'
#' @param chromosomes.as.numbers If TRUE, replace chromosome names with numbers for readability (default: FALSE)
#'
#' @param show.chromosome.names If TRUE, display chromosome names on the x axis (default: TRUE)
#'
#' @param x.axis.title Title to display on the x axis (default: NULL, i.e. no title)
#'
#' @param default.point.color Default color for a track when not specified in track data,
#' either a string (e.g. "grey20") or a vector of alternating colors (e.g. c("red", "blue") for two colors) (default: c("dodgerblue3", "darkgoldenrod2"))
#'
#' @param default.bg.color Default background color when not specified in track data,
#' either a string (e.g. "grey20") or a vector of alternating colors (e.g. c("grey80", "white") for two colors) (default: c("grey85", "white"))
#'
#' @param default.point.size Default point size for a track when not specified in track data (default: 0.5)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @examples
#' # Manhattan plot showing female-specific and male-specific SNPs on two tracks
#' plot_manhattan("psass_window.tsv",
#'                tracks = list(manhattan_track("Snps_females", point.color = c("firebrick1", "firebrick3")),
#'                              manhattan_track("Snps_males", point.color = c("dodgerblue1", "dodgerblue3"))),
#'                output.file = "manhattan.png",
#'                chromosomes.as.numbers = TRUE, show.chromosome.names = TRUE, x.axis.title = "Chromosome")
#'

plot_manhattan <- function(input.file.path, tracks,
                           chromosomes.file.path = NULL, detect.chromosomes = TRUE, unplaced.label = "Unplaced",
                           output.file = NULL, width = 14, track.height = 6, res = 300,
                           chromosomes.as.numbers = FALSE, show.chromosome.names = TRUE,
                           x.axis.title = NULL,
                           default.point.color = c("dodgerblue3", "darkgoldenrod2"),
                           default.bg.color = c("grey85", "white"),
                           default.point.size = 0.5, default.ylim = NULL) {


    # Load chromosome names (return NULL of no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes.file.path)

    # Load genomic data
    data <- load_genome_input(input.file.path, chromosomes = chromosomes, detect.chromosomes = detect.chromosomes, unplaced.label = unplaced.label)

    manhattan_plot <- draw_manhattan_plot(data$data,
                                          data$lengths,
                                          tracks,
                                          output.file = output.file,
                                          width = width,
                                          track.height = track.height,
                                          res = res,
                                          chromosomes.as.numbers = chromosomes.as.numbers,
                                          show.chromosome.names = show.chromosome.names,
                                          x.axis.title = x.axis.title,
                                          default.point.color = default.point.color,
                                          default.bg.color = default.bg.color,
                                          default.point.size = default.point.size,
                                          default.ylim = default.ylim)
}
