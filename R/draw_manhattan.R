#' @title Draw manhattan plot
#'
#' @description Generate a manhattan plot with multiple tracks for the entire genome
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param contig.lengths Contig lengths vector (e.g. from the output of the \code{\link{load_genome_input}} function)
#'
#' @param tracks List of tracks to plot. Tracks can be generated with the \code{\link{manhattan_track}} function
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
#' @return Combined plot data (ggplot object)
#'
#' @examples
#' # Manhattan plot showing female-specific and male-specific SNPs on two tracks
#' genomic_data <- load_genome_input("psass_window.tsv")
#' draw_manhattan_plot(genomic_data$data, genomic_data$lengths,
#'                     tracks = list(manhattan_track("Snps_females", point.color = c("firebrick1", "firebrick3")),
#'                                   manhattan_track("Snps_males", point.color = c("dodgerblue1", "dodgerblue3"))),
#'                     output.file = "manhattan.png",
#'                     chromosomes.as.numbers = TRUE, show.chromosome.names = TRUE, x.axis.title = "Chromosome")
#'

draw_manhattan_plot <- function(data, contig.lengths, tracks,
                                output.file = NULL, width = 14, track.height = 6, res = 300,
                                chromosomes.as.numbers = FALSE, show.chromosome.names = TRUE,
                                x.axis.title = NULL,
                                default.point.color = c("dodgerblue3", "darkgoldenrod2"),
                                default.bg.color = c("grey85", "white"),
                                default.point.size = 0.5, default.ylim = NULL) {

    # Compute increment to add to each contig (i.e. cumulative length of each contig before this one)
    increments <- setNames(c(0, cumsum(head(contig.lengths, -1))), names(contig.lengths))

    # Adjust x-axis position for each point in the data based on increments
    data$Position_plot = data$Position_plot + increments[data$Contig_plot]

    # Create data for background rectangles (for alternating background color)
    backgrounds <- data.frame(contig = names(increments), start = increments, end = increments + contig.lengths)

    # Initialize list of plots
    n_tracks <- length(tracks)
    plots <- rep(list(NULL), n_tracks)

    x.labels.angle = 90  # Vertical labels when full chromosome names
    if (chromosomes.as.numbers) {x.labels.angle = 0}  # Horizontal labels when numbers

    # Draw specified tracks
    bottom_track <- FALSE
    for (i in c(1:n_tracks)) {

        if (i == n_tracks) bottom_track <- TRUE  # For x-axis labels and title

        # Assign default values to track properties if not specified by user
        tracks[[i]] <- assign_manhattan_track_default(tracks[[i]], default.point.color = default.point.color, default.bg.color = default.bg.color,
                                                      default.point.size = default.point.size, default.ylim = default.ylim)

        # Generate track data
        track_data <- create_manhattan_track_data(data, tracks[[i]], chromosomes.as.numbers = chromosomes.as.numbers)

        # Generate background data
        track_background_data <- create_manhattan_background_data(backgrounds, tracks[[i]])

        # Generate track plot
        plots[[i]] <- plot_track_manhattan(track_data, track_background_data, tracks[[i]], bottom.track = bottom_track,
                                           show.chromosome.names = show.chromosome.names, x.labels.angle = x.labels.angle,
                                           x.axis.title = x.axis.title)

    }

    # Combine all tracks in a single plot
    combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")

    # Output to file if specified or print in current R device otherwise
    if (!is.null(output.file)) {

        ggplot2::ggsave(output.file, plot = combined, width = width, height = track.height * n_tracks, dpi = res)

    } else {

        print(combined)

    }

    return(combined)
}



#' @title Create manhattan track object
#'
#' @description Generate an object storing all properties for a manhattan track
#'
#' @param metric Metric represented by the track, which should correspond to a column name in the data frame used as input
#' data in the plot (output of the \code{\link{load_genome_input}} function)
#'
#' @param label Track label, NULL to set the label to the metric name (default: NULL)
#'
#' @param point.color Point color, either a string (e.g. "grey20") or a vector of alternating colors (e.g. c("red", "blue") for two colors)
#'
#' @param bg.color Background color, either a string (e.g. "white") or a vector of alternating colors (e.g. c("grey85", "white") for two colors)
#'
#' @param point.size Point size, directly passed to ggplot
#'
#' @param ylim Vector of y-axis limits for the track; if NULL, infer directly from data
#'
#' @return A named list with the value of each track property
#'
#' @examples
#' # Single metric
#' track_data <- manhattan_track("Fst", color = c("blue", "yellow"), point.size = 0.75)

manhattan_track <- function(metric, label = NULL, point.color = NULL, bg.color = NULL, point.size = NULL, ylim = NULL) {

    # Set a default label if not specified
    if (is.null(label)) { label <- metric }

    # Create track object
    track_info <- list(metric=metric, label=label, point.color=point.color, bg.color=bg.color, point.size=point.size, ylim=ylim)

    return(track_info)
}



#' @title Assign default values to a manhattan track object
#'
#' @description Assign default values to all properties for which the value was not specified by the user (e.g. value is NULL)
#'
#' @param track A track object generated with the \code{\link{manhattan_track}} function)
#'
#' @param default.point.color Default point color when not specified in track data (default: c("dodgerblue3", "darkgoldenrod2"))
#'
#' @param default.bg.color Default background color when not specified in track data (default: c("white", "grey85"))
#'
#' @param default.point.size Default point size for a track when not specified in track data (default: 0.5)
#'
#' @param default.ylim Default y-axis limits for a track when not specified in track data (default: NULL, i.e. infer from data)
#'
#' @return A track object with default values for properties not specified by the user
#'
#' @examples
#' track_data <- assign_manhattan_track_default(track_data, default.point.size = 0.75)
#'

assign_manhattan_track_default <- function(track, default.point.color = c("dodgerblue3", "darkgoldenrod2"),
                                           default.bg.color = c("white", "grey85"),
                                           default.point.size = 0.5, default.ylim = NULL) {

    if (is.null(track$point.color)) { track$point.color <- default.point.color }
    if (is.null(track$bg.color)) { track$bg.color <- default.bg.color }
    if (is.null(track$point.size)) { track$point.size <- default.point.size }
    if (is.null(track$ylim)) { track$ylim <- default.ylim }

    return(track)

}



#' @title Create manhattan track data
#'
#' @description Create input data frame for the \code{\link{plot_track_manhattan}} function from the genomic data and the track information
#'
#' @param data Genomic data (e.g. result of PSASS or RADSex loaded with the \code{\link{load_genome_input}} function)
#'
#' @param track Track object for the current plot, generated with the \code{\link{manhattan_track}} function
#'
#' @param chromosomes.as.numbers If TRUE, replace chromosome names with numbers for readability (default: FALSE)
#'
#' @return A data frame with columns:
#' Contig_plot | Position_plot | Metric | Point Color
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' track <- manhattan_track("Fst")
#'
#' track_data <- create_manhattan_track_data(genomic_data, track)
#'

create_manhattan_track_data <- function(data, track, chromosomes.as.numbers = FALSE) {

    # Extract required columns and create color columns
    track_data <- data[, c("Contig_plot", "Position_plot", track$metric)]
    names(track_data)[3] <- "Value"
    # Sort data by Contig then Position
    track_data <- track_data[order(track_data$Contig_plot, track_data$Position_plot), ]
    # Assign point color
    contigs <- unique(track_data$Contig_plot)
    n_contigs <- length(contigs)
    points_palette <- setNames(rep(track$point.color, n_contigs)[1:n_contigs], contigs)
    track_data$Color <- points_palette[track_data$Contig_plot]

    if (chromosomes.as.numbers) {
        chromosome_numbers <- setNames(c(seq(1, n_contigs - 1), 'U'), contigs)
        track_data$Contig_plot <- chromosome_numbers[track_data$Contig_plot]
    }

    return(track_data)
}



#' @title Create manhattan background data
#'
#' @description Create background data frame for the \code{\link{plot_track_manhattan}} function
#'
#' @param data Data frame with contig, start, and end for each background rectangle
#'
#' @param track Track object for the current plot, generated with the \code{\link{manhattan_track}} function
#'
#' @return A data frame with columns:
#' contig | start | end | color
#'
#' @examples
#' track <- manhattan_track("Fst")
#' background_data <- create_manhattan_background_data(background_data, track)
#'

create_manhattan_background_data <- function(data, track) {

    background_data <- data

    # Assign point and background color
    n_contigs <- nrow(data)
    bg_palette <- setNames(rep(track$bg.color, n_contigs)[1:n_contigs], data$contig)
    background_data$color <- bg_palette[background_data$contig]

    return(background_data)
}
