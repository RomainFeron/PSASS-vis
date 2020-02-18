#' @title Plot region
#'
#' @description Generate a linear plot with multiple tracks for a specified genomic region from a genomic data file
#'
#' @param input.file.path Path to a genomic data input file (e.g. result of PSASS or RADSex)
#'
#' @param region Region to plot, with syntax "Contig" or "Contig:start-end"
#'
#' @param tracks List of tracks to plot. Tracks can be generated with the \code{\link{region_track}} function
#'
#' @param chromosomes.file.path Path to a tabulated file specifying chromosome names (default: NULL)
#'
#' @param detect.chromosomes If TRUE, will consider contigs starting with LG, CH, or NC as chromosomes if no chromosomes were specified (default: TRUE)
#'
#' @param output.file Path to an output file for the generated region plot, or NULL to plot in the current R device (default: NULL)
#'
#' @param width Plot width when plotting to an output file, in inches (default: 12)
#'
#' @param track.height Height of a single track when plotting to an output file, in inches (default: 4)
#'
#' @param res Image resolution when plotting to an output file, in dpi (default: 300)
#'
#' @param default.color Default color for a track when not specified in track data (default: "grey20")
#'
#' @param default.alpha Default alpha value for a track when not specified in track data (default: 1)
#'
#' @param default.type Default plot type for a track when not specified in track data (default: "ribbon")
#'
#' @param default.point.size Default point size for a track of type "points" when not specified in track data (default: 0.5)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @param default.major.lines.y If TRUE, reference lines will be plotted for the y axis if not specified in track data (default: TRUE)
#'
#' @param default.major.lines.x If TRUE, reference lines will be plotted for the x axis if not specified in track data (default: FALSE)
#'
#' @param default.legend.position Default legend position when not specified in track data (default: "none")
#'
#' @return Combined plot data (ggplot object)
#'
#' @examples
#'
#' # Plotting an FST track and a combined pool-specific SNPs track
#' # for an entire chromosome to the default R device
#'
#' plot_region("data/psass_window.tsv", "Chr01",
#'             tracks = list(region_track("Fst", label = expression("F"["ST"])),
#'                           region_track(c("Snps_females", "Snps_males"), label = "Pool-specific SNPs", color = c("firebrick2", "dodgerblue3"), alpha=0.6)),
#'             chromosomes.file.path = "data/chromosomes.tsv")
#'

plot_region <- function(input.file.path, region, tracks,
                        chromosomes.file.path = NULL, detect.chromosomes = TRUE,
                        output.file = NULL, width = 12, track.height = 4, res = 300,
                        default.color = "grey20", default.alpha = 1, default.type = "ribbon", default.point.size = 0.5, default.ylim = NULL,
                        default.major.lines.y = TRUE, default.major.lines.x = FALSE, default.legend.position = "none") {


    # Load chromosome names (return NULL of no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes.file.path)

    # Load genomic data
    data <- load_genome_input(input.file.path, chromosomes = chromosomes, detect.chromosomes = detect.chromosomes, unplaced.label = "Unplaced")

    # Plot genomic region
    region_plot <- draw_region(data$data,
                               data$lengths,
                               region,
                               tracks,
                               output.file = output.file,
                               width = width,
                               track.height = track.height,
                               res = res,
                               default.color = default.color,
                               default.alpha = default.alpha,
                               default.type = default.type,
                               default.point.size = default.point.size,
                               default.ylim = default.ylim,
                               default.major.lines.y = default.major.lines.y,
                               default.major.lines.x = default.major.lines.x,
                               default.legend.position = default.legend.position)
}
