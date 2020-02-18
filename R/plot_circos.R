#' @title Plot circos
#'
#' @description Generate a circos plot with multiple tracks from a genomic data file
#'
#' @param input.file.path Path to a genomic data input file (e.g. result of PSASS or RADSex)
#'
#' @param tracks List of tracks to plot. Tracks can be generated with the \code{\link{region_track}} function
#'
#' @param chromosomes.file.path Path to a tabulated file specifying chromosome names (default: NULL)
#'
#' @param detect.chromosomes If TRUE, will consider contigs starting with LG, CH, or NC as chromosomes if no chromosomes were specified (default: TRUE)
#'
#' @param unplaced.label Label for unplaced contigs (default: "Unplaced")
#'
#' @param output.file Path to an output file for the generated circos plot, or NULL to plot in the current R device (default: NULL)
#'
#' @param width Plot width when plotting to an output file, in pixel (default: 2400)
#'
#' @param height Plot height when plotting to an output file, in pixel (default: 2400)
#'
#' @param res Image resolution when plotting to an output file, in ppi (default: 120)
#'
#' @param default.color Default color for a track when not specified in track data (default: "grey20")
#'
#' @param default.point.size Default point size for a track when not specified in track data (default: 0.25)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @param default.bg.color Default background color for standard sectors in the track (default: "white")
#'
#' @param default.bg.highlight.color Default background color for highlighted sectors in the track (default: "grey80")
#'
#' @param highlight Vector of names of sectors to highlight
#'
#' @param sector.titles.expand Manually set the space between sector titles and x-axis as a multiple of ymax (default: NULL)
#'
#' @examples
#' plot_circos("data/psass_window.tsv",
#'             tracks = list(track("Fst", label = expression("F"["ST"])),
#'                           track(c("Snps_females", "Snps_males"), label = "Pool-specific SNPs", color = c("firebrick2", "dodgerblue3")),
#'                           track(c("Depth_ratio"), label = "Depth ratio", color = "grey50")),
#'             output.file = "circos.png")
#'

plot_circos <- function(input.file.path, tracks,
                        chromosomes.file.path = NULL, detect.chromosomes = TRUE, unplaced.label = "Unplaced",
                        highlight = NULL,
                        output.file = NULL, width = 2400, height = 2400, res = 120,
                        default.color = "grey20", default.point.size = 0.25, default.ylim = NULL,
                        default.bg.color = "white", default.bg.highlight.color = "grey80",
                        sector.titles.expand = NULL) {

    # Load chromosome names (return NULL of no chromosomes file)
    chromosomes <- load_chromosome_names(chromosomes.file.path)

    # Load genomic data
    data <- load_genome_input(input.file.path, chromosomes = chromosomes, detect.chromosomes = detect.chromosomes, unplaced.label = unplaced.label)

    # Draw the plot
    draw_circos(data$data,
                data$lengths,
                tracks,
                highlight = highlight,
                output.file = output.file,
                width = width,
                height = height,
                res = res,
                default.color = default.color,
                default.point.size = default.point.size,
                default.ylim = default.ylim,
                default.bg.color = default.bg.color,
                default.bg.highlight.color = default.bg.highlight.color,
                sector.titles.expand = sector.titles.expand
    )
}
