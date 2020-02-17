#' @title Draw circos plot
#'
#' @description Generate a circular plot with multiple track for the entire genome
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param contig.lengths Contig lengths from the output of the \code{\link{load_genome_input}} function
#'
#' @param tracks List of tracks to plot. Tracks can be generated with the \code{\link{region_track}} function
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
#' genomic_data <- load_genome_input("psass_window.tsv")
#' draw_circos(genomic_data$data, genomic_data$lengths,
#'             tracks = list(track("Fst", label = expression("F"["ST"])),
#'                           track(c("Snps_females", "Snps_males"), label = "Pool-specific SNPs", color = c("firebrick2", "dodgerblue3")),
#'                           track(c("Depth_ratio"), label = "Depth ratio", color = "grey50")),
#'             output.file = "circos.png")
#'

draw_circos <- function(data, contig.lengths, tracks,
                        highlight = NULL,
                        output.file = NULL, width = 2400, height = 2400, res = 120,
                        default.color = "grey20", default.point.size = 0.25, default.ylim = NULL,
                        default.bg.color = "white", default.bg.highlight.color = "grey80",
                        sector.titles.expand = NULL) {

    # Check that sectors to highlight exist
    highlight <- check_highlight_sectors(data, contig.lengths, highlight)

    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, length(contig.lengths)), contig.lengths), ncol=2)

    # Setup gaps between sectors
    gaps <- rep(1, nrow(sector_lengths))
    gaps[length(gaps)] <- 12  # Bigger gap after the last sector to leave some space for track titles

    # Set a small offset to starting angle to make room for axes labels
    starting_angle_offset <- 7.5

    # Calculate width of track and sector.title.expand factor based on number of tracks
    n_tracks <- length(tracks)
    track_height <- 0.25  # Width of 0.25 if single track is plotted
    if (n_tracks > 1) { track_height <- 0.5 / n_tracks }  # For multiple tracks, total width 0.5 (half the circo's width)
    if (is.null(sector.titles.expand)) sector.titles.expand <- 1.1 + n_tracks / 10

    # Reset circos parameters
    circlize::circos.clear()

    # Setup circos parameters
    circlize::circos.par("track.height" = track_height,
                         "start.degree" = 90 - starting_angle_offset,
                         "gap.degree" = gaps,
                         "cell.padding" = c(0.01, 0.1, 0.01, 0.1),
                         "points.overflow.warning" = FALSE,
                         "track.margin" = c(0, 0.015),
                         "unit.circle.segments" = 100)

    # Open output file if specified
    if (!is.null(output.file)) { png(output.file, width=width, height=height,  res=res) }

    # Initialize circos plot
    circlize::circos.initialize(factors = names(contig.lengths), xlim = sector_lengths, sector.width = contig.lengths)

    # Draw specified tracks
    top_track <- TRUE
    for (i in c(1:length(tracks))) {

        # Assign default values to track properties when not set by user
        tracks[[i]] <- assign_circos_track_default(tracks[[i]], default.color = default.color, default.point.size = default.point.size,
                                                   default.ylim = default.ylim, default.bg.color = default.bg.color,
                                                   default.bg.highlight.color = default.bg.highlight.color)

        # Setup sector background colors for the track
        if (length(tracks[[i]]$bg.color) == 1) {  # Background colors defined as singe values for "bg.color" and "highlight.bg.color"

            tracks[[i]]$bg_colors <- rep(tracks[[i]]$bg.color, length(contig.lengths))
            tracks[[i]]$bg_colors[which(names(contig.lengths) %in% highlight)] <- tracks[[i]]$bg.highlight.color

        } else {  # Background colors defined as vector of values for "bg.color"

            tracks[[i]]$bg_colors <- tracks[[i]]$bg.color

        }

        # Create input data for track
        track_data <- create_circos_track_data(data, tracks[[i]])

        # Plot a single track
        plot_circos_track(track_data, tracks[[i]],
                          top.track = top_track, sector.titles.expand = sector.titles.expand,
                          first_sector = names(contig.lengths)[1])

        top_track <- FALSE
    }

    # Close output file
    if (!is.null(output.file)) { dev.off() }
}



#' @title Check highlight sectors
#'
#' @description Check that sectors to be highlighted exist, allowing both contig names and chromosomes names
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param contig.lengths Contig lengths from the output of the \code{\link{load_genome_input}} function
#'
#' @param highlight Vector of names of sectors to highlight
#'
#' @return A clean vector of names for sectors to highlight
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#'
#' highlight <- check_highlight_sectors(genomic_data$data, genomic_data$lengths, c("Chr01", "NC_002364.1"))
#'

check_highlight_sectors <- function(data, contig.lengths, highlight) {

    if (is.null(highlight)) {return(c())}

    for (i in 1:length(highlight)) {

        if (!(highlight[i] %in% names(contig.lengths))) {  # Sector to highlight not found in list of sector lengths

            contigs <- setNames(unique(data$Contig_plot), unique(data$Contig))

            if (highlight[i] %in% names(contigs)) {  # Look for specified sector to highlight in original contig names

                highlight[i] <- contigs[highlight[i]]

            } else {

                print(paste0("Could not find sector to highlight \"", highlight[i], "\"."))
                exit(1)

            }
        }
    }

    return(highlight)
}



#' @title Create circos track data
#'
#' @description Create input data frame for the \code{\link{plot_track_circos}} function from the genomic data and the track information
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param track Track object for the current plot, generated with the \code{\link{circos_track}} function
#'
#' @return A data frame with columns:
#' Contig_plot | Position_plot | Metric 1 | Metric 1 colors | ... | Metric N | Metric N colors
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' track <- circos_track("Fst")
#'
#' region_data <- create_circos_track_data(genomic_data, track)
#'

create_circos_track_data <- function(data, track) {

    # Extract required columns and create color columns
    track_data <- data[, c("Contig_plot", "Position_plot", track$metrics)]
    # Combine data for multiple metrics
    track_data <- reshape2::melt(track_data, measure.vars = track$metrics, variable.name="Color")
    # Assign color to data points
    track_data$Color <- setNames(track$color, track$metrics)[track_data$Color]
    # Sort data by Contig then Position
    track_data <- track_data[order(track_data$Contig_plot, track_data$Position_plot), ]

    return(track_data)
}



#' @title Create circos track object
#'
#' @description Generate an object storing all properties for a circos track
#'
#' @param metrics Metrics included in the track. Metrics should correspond to column names in the data frame used as input
#' data in the plot (output of the \code{\link{load_genome_input}} function). Possible values: a string if the track includes
#' a single metric (e.g. "Fst"), or a vector if the track includes several metrics (e.g. c("Snps_females", "Snps_males"))
#'
#' @param label Track label, NULL to set the label to the metric name for single-metric tracks. Label has to be specified for multi-metrics tracks (default: NULL)
#'
#' @param color Track color. Values can be a string (e.g. "grey20") and will then be applied to all metrics, or a vector of size equal to the number of metrics
#' (e.g. c("red", "blue") for two metrics)
#'
#' @param point.size Point size for plots of type "points". Values can be a float (e.g. 0.5) and will then be applied to all metrics,
#' or a vector of size equal to the number of metrics (e.g. c(1, 1.5, 3) for three metrics)
#'
#' @param ylim Vector of y-axis limits for the track; if NULL, infer directly from data
#'
#' @param bg.color Background color for standard sectors in the track (single value or vector of length equal to the number of sectors in the plot)
#'
#' @param bg.highlight.color Background color for highlighted sectors in the track (single value)
#'
#' @return A named list with the value of each track property
#'
#' @examples
#' # Single metric
#' track_data <- circos_track("Fst", color = "grey70", point.size = 0.75)
#'
#' # Multiple metrics
#' track_data <- circos_track(c("Snp_females", "Snp_males"), color = c("firebrick2", "dodgerblue3"))
#'

circos_track <- function(metrics, label = NULL, color = NULL, point.size = NULL, ylim = NULL,
                         bg.color = NULL, bg.highlight.color = NULL) {

    n_metrics <- length(metrics)

    if (is.null(label)) {

        if (n_metrics == 1) {

            label <- metrics  # Set a default label if not specified

        } else {

            stop("Label required for multi-metrics track")

        }
    }

    # Handle multi-values properties for multiple metrics (i.e. metric-specific properties)
    if (n_metrics > 1) {

        # For each option, assign value to each metric if single value was defined in multi-metrics track
        if (length(color) == 1) { color <- rep(color, n_metrics) }
        if (length(point.size) == 1) { point.size <- rep(point.size, n_metrics) }

    }

    # Single-value properties (not metric-specific)
    ylim <- ylim
    bg.color <- bg.color
    bg.highlight.color <- bg.highlight.color

    track_info <- list(metrics=metrics, label=label, color=color, point.size=point.size, ylim = ylim,
                       bg.color = bg.color, bg.highlight.color = bg.highlight.color)

    return(track_info)
}



#' @title Assign default values to a circos track object
#'
#' @description Assign default values to all properties for which the value was not specified by the user (e.g. value is NULL)
#'
#' @param track A track object generated with the \code{\link{circos_track}} function)
#'
#' @param default.color Default color for a track when not specified in track data (default: "grey20")
#'
#' @param default.point.size Default point size for a track of type "points" when not specified in track data (default: 0.5)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @param default.bg.color Default background color for standard sectors in the track (default: "white")
#'
#' @param default.bg.highlight.color Default background color for highlighted sectors in the track (default: "grey80")
#'
#' @return A track object with default values for properties not specified by the user
#'
#' @examples
#' track_data <- assign_circos_track_default(track_data, default.alpha = 0.75)
#'

assign_circos_track_default <- function(track, default.color = "grey20", default.point.size = 0.5, default.ylim = NULL,
                                        default.bg.color = "white", default.bg.highlight.color = "grey80") {

    n_metrics <- length(track$metrics)

    # Metrics-specific options (create vector)
    if (is.null(track$color)) { track$color <- rep(default.color, n_metrics) }
    if (is.null(track$point.size)) { track$point.size <- rep(default.point.size, n_metrics) }

    # Track-specific options (single value)
    if (is.null(track$ylim)) { track$ylim <- default.ylim }
    if (is.null(track$bg.color)) { track$bg.color <- default.bg.color }
    if (is.null(track$bg.highlight.color)) { track$bg.highlight.color <- default.bg.highlight.color }

    return(track)

}

