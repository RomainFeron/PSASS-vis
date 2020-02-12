
draw_region <- function(data, contig_lengths, region,
                        tracks = NULL,
                        default.color = "grey20", default.alpha = 1, default.type = "ribbon", default.point.size = 0.5, default.ylim=NULL,
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

        if (i == n_tracks) bottom_track <- TRUE

        tracks[[i]] <- assign_track_default(tracks[[i]], default.color = default.color,
                                            default.alpha = default.alpha, default.type = default.type, default.point.size = default.point.size,
                                            default.major.lines.y = default.major.lines.y, default.major.lines.x = default.major.lines.x,
                                            default.legend.position = default.legend.position, default.ylim = default.ylim)

        track_data <- create_region_track_data(data, region_info, tracks[[i]])

        plots[[i]] <- track_region(track_data, region_info, tracks[[i]], bottom.track = bottom_track)

    }

    combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")

    print(combined)

    if (!is.null(output.file)) {

        ggplot2::ggsave(output.file, plot = combined, width = width, height = track.height * n_tracks, dpi = res)

    } else {

        print(combined)

    }

    return(combined)
}



create_region_track_data <- function(data, region_info, track) {

    # Extract region from data
    if (region_info[[1]] %in% unique(data$Contig)) {

        data = subset(data, data$Contig == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    } else if (region_info[[1]] %in% unique(data$Contig_plot)) {

        data = subset(data, data$Contig_plot == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    }

    # Extract required columns and create color columns
    track_data <- data[, c("Position")]
    for (i in 1:length(track$metrics)) { track_data = cbind(track_data, data[, c(track$metrics[i])], rep(track$color[i], nrow(data))) }

    return(track_data)
}



track <- function(metrics, label = NULL, color = NULL, alpha = NULL, type = NULL, point.size = NULL,
                  major.lines.y = NULL, major.lines.x = NULL, legend.position = NULL, ylim=NULL) {

    n_metrics <- length(metrics)

    if (is.null(label)) {

        if (n_metrics ==1) {

            label <- metrics  # Set a default label if not specified

        } else {

            stop("Label required for multi-metrics track")

        }
    }

    if (n_metrics > 1) {

        # For each option, assign value to each metric if single value was defined in multi-metrics track
        if (length(color) == 1) { color <- rep(color, n_metrics) }
        if (length(alpha) == 1) { alpha <- rep(alpha, n_metrics) }
        if (length(type) == 1) { type <- rep(type, n_metrics) }
        if (length(point.size) == 1) { point.size <- rep(point.size, n_metrics) }
        if (length(major.lines.y) == 1) { major.lines.y <- rep(major.lines.y, n_metrics) }
        if (length(major.lines.x) == 1) { major.lines.x <- rep(major.lines.x, n_metrics) }
        if (length(legend.position) == 1) { legend.position <- rep(legend.position, n_metrics) }
        if (length(ylim) == 1) { ylim <- rep(ylim, n_metrics) }

    }

    track_info <- list(metrics=metrics, label=label, color=color, alpha=alpha, type=type, point.size=point.size,
                       major.lines.y = major.lines.y, major.lines.x = major.lines.x, legend.position = legend.position,
                       ylim = ylim)

    return(track_info)
}


assign_track_default <- function(track, default.color = "grey20", default.alpha = 1, default.type = "ribbon",
                                 default.point.size = 0.5, default.ylim=NULL,
                                 default.major.lines.y = TRUE, default.major.lines.x = FALSE,
                                 default.legend.position = "right") {


    n_metrics <- length(track$metrics)

    # Metrics-specific options (create vector)
    if (is.null(track$color)) { track$color = rep(default.color, n_metrics) }
    if (is.null(track$alpha)) { track$alpha = rep(default.alpha, n_metrics) }
    if (is.null(track$type)) { track$type = rep(default.type, n_metrics) }
    if (is.null(track$point.size)) { track$point.size = rep(default.point.size, n_metrics) }
    if (is.null(track$legend.position)) { track$legend.position = rep(default.legend.position, n_metrics) }
    if (is.null(track$ylim)) { track$ylim = rep(default.ylim, n_metrics) }

    # Track-specific options (single value)
    if (is.null(track$major.lines.y)) { track$major.lines.y = default.major.lines.y }
    if (is.null(track$major.lines.x)) { track$major.lines.x = default.major.lines.x }

    return(track)

}