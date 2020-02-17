#' @title Draw region plot
#'
#' @description Generate a linear plot with multiple tracks for a specified genomic region
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param contig_lengths Contig lengths from the output of the \code{\link{load_genome_input}} function
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
#' # for the first 3Mb of chromosome 1 to the default R device
#'
#' genomic_data <- load_genome_input("psass_window.tsv")
#'
#' draw_region(genomic_data$data, genomic_data$lengths, "Chr01:0-3000000",
#'             tracks = list(region_track("Fst", label = expression("F"["ST"])),
#'                           region_track(c("Snps_females", "Snps_males"), label = "Pool-specific SNPs", color = c("firebrick2", "dodgerblue3"), alpha=0.6)))
#'

draw_region <- function(data, contig_lengths, region,
                        tracks = NULL,
                        default.color = "grey20", default.alpha = 1, default.type = "ribbon", default.point.size = 0.5, default.ylim = NULL,
                        default.major.lines.y = TRUE, default.major.lines.x = FALSE, default.legend.position = "right",
                        output.file = NULL, width = 12, track.height = 4, res = 300) {

    # Add original contig names to contig lengths so that user can specify both chromosome names or contig names in region
    contig_lengths <- c(contig_lengths, setNames(unique(data$Length), unique(data$Contig)))

    # Get contig, min, and max from the region string
    region_info <- parse_region(region, contig_lengths)

    # Initialize list of plots
    n_tracks <- length(tracks)
    plots <- rep(list(NULL), n_tracks)

    # Draw specified tracks
    bottom_track <- FALSE
    for (i in c(1:n_tracks)) {

        if (i == n_tracks) bottom_track <- TRUE  # For x-axis labels and title

        # Assign default values to track properties if not specified by user
        tracks[[i]] <- assign_region_track_default(tracks[[i]], default.color = default.color,
                                                   default.alpha = default.alpha, default.type = default.type, default.point.size = default.point.size,
                                                   default.major.lines.y = default.major.lines.y, default.major.lines.x = default.major.lines.x,
                                                   default.legend.position = default.legend.position, default.ylim = default.ylim)

        # Generate track data
        track_data <- create_region_track_data(data, region_info, tracks[[i]])

        # Generate track plot
        plots[[i]] <- plot_track_region(track_data, region_info, tracks[[i]], bottom.track = bottom_track)

    }

    # Combine all tracks in a single plot
    combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")

    # Output to file if specified or print in current R device otherwise
    if (!is.null(output.file)) {

        ggplot2::ggsave(output.file, plot = combined, width = width, height = track.height * n_tracks, dpi = res)

    } else {

        print(combined)

    }

    # Return combined ggplot object
    return(combined)
}



#' @title Create region track data
#'
#' @description Create input data frame for the \code{\link{plot_track_region}} function from the genomic data and the track information
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param region.info Information on the region to plot, output of the \code{\link{parse_region}} function
#'
#' @param track Track object for the current plot, generated with the \code{\link{region_track}} function
#'
#' @return A data frame with columns:
#' Position | Metric 1 | Metric 1 colors | ... | Metric N | Metric N colors
#'
#' @examples
#'
#' genomic_data <- load_genome_input("psass_window.tsv")
#' region_info <- parse_region("Chr01:0-1500000")
#' track <- region_track("Fst")
#'
#' region_data <- create_region_track_data(genomic_data, region_info, track)
#'

create_region_track_data <- function(data, region.info, track) {

    # Extract data points for the region from genomic data
    if (region_info[[1]] %in% unique(data$Contig)) {

        data <- subset(data, data$Contig == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    } else if (region_info[[1]] %in% unique(data$Contig_plot)) {

        data <- subset(data, data$Contig_plot == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    }

    # Extract required columns and create color columns
    track_data <- data[, c("Position")]
    for (i in 1:length(track$metrics)) { track_data = cbind(track_data, data[, c(track$metrics[i])], rep(track$color[i], nrow(data))) }

    return(track_data)
}



#' @title Create region track object
#'
#' @description Generate an object storing all properties for a region track
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
#' @param alpha Track alpha value. Values can be a float (e.g. 0.5) and will then be applied to all metrics, or a vector of size equal to the number of metrics
#' (e.g. c(1, 0.5, 0.75) for three metrics)
#'
#' @param type Track plot type. Possible values: "ribbon" or "points". Values can be a string (e.g. "ribbon") and will then be applied to all metrics,
#' or a vector of size equal to the number of metrics (e.g. c("points", "ribbon") for two metrics)
#'
#' @param point.size Point size for plots of type "points". Values can be a float (e.g. 0.5) and will then be applied to all metrics,
#' or a vector of size equal to the number of metrics (e.g. c(1, 1.5, 3) for three metrics)
#'
#' @param major.lines.y Plot reference lines for the y axis, directly passed to "major.lines.y" in \code{\link{ggplot2::theme}} (single value)
#'
#' @param major.lines.x Plot reference lines for the x axis, directly passed to "major.lines.x" in \code{\link{ggplot2::theme}} (single value)
#'
#' @param legend.position Position of the legend for the track, directly passed to "legend.position" in \code{\link{ggplot2::theme}} (single value)
#'
#' @param ylim Vector of y-axis limits for the track; if NULL, infer directly from data
#'
#' @return A named list with the value of each track property
#'
#' @examples
#'
#' # Single metric
#' track_data <- region_track("Fst", color = "grey70", type = "points", point.size = 0.75)
#'
#' # Multiple metrics
#' track_data <- region_track(c("Snp_females", "Snp_males"), color = c("firebrick2", "dodgerblue3"), alpha = 0.5)
#'

region_track <- function(metrics, label = NULL, color = NULL, alpha = NULL, type = NULL, point.size = NULL,
                         major.lines.y = NULL, major.lines.x = NULL, legend.position = NULL, ylim = NULL) {

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
        if (length(alpha) == 1) { alpha <- rep(alpha, n_metrics) }
        if (length(type) == 1) { type <- rep(type, n_metrics) }
        if (length(point.size) == 1) { point.size <- rep(point.size, n_metrics) }

    }

    # Single-value properties (not metric-specific)
    major.lines.y <- major.lines.y
    major.lines.x <- major.lines.x
    legend.position <- legend.position
    ylim <- ylim

    track_info <- list(metrics=metrics, label=label, color=color, alpha=alpha, type=type, point.size=point.size,
                       major.lines.y = major.lines.y, major.lines.x = major.lines.x, legend.position = legend.position,
                       ylim = ylim)

    return(track_info)
}



#' @title Assign default values to a region track object
#'
#' @description Assign default values to all properties for which the value was not specified by the user (e.g. value is NULL)
#'
#' @param track A track object generated with the \code{\link{region_track}} function)
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
#' @return A track object with default values for properties not specified by the user
#'
#' @examples
#'
#' track_data <- assign_track_default(track_data, default.alpha = 0.75)
#'

assign_region_track_default <- function(track, default.color = "grey20", default.alpha = 1, default.type = "ribbon",
                                        default.point.size = 0.5, default.ylim = NULL,
                                        default.major.lines.y = TRUE, default.major.lines.x = FALSE,
                                        default.legend.position = "right") {

    n_metrics <- length(track$metrics)

    # Metrics-specific options (create vector)
    if (is.null(track$color)) { track$color <- rep(default.color, n_metrics) }
    if (is.null(track$alpha)) { track$alpha <- rep(default.alpha, n_metrics) }
    if (is.null(track$type)) { track$type <- rep(default.type, n_metrics) }
    if (is.null(track$point.size)) { track$point.size <- rep(default.point.size, n_metrics) }

    # Track-specific options (single value)
    if (is.null(track$legend.position)) { track$legend.position <- default.legend.position }
    if (is.null(track$major.lines.y)) { track$major.lines.y <- default.major.lines.y }
    if (is.null(track$major.lines.x)) { track$major.lines.x <- default.major.lines.x }
    if (is.null(track$ylim)) { track$ylim <- default.ylim }

    return(track)

}