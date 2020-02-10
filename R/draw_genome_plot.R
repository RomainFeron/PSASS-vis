
draw_circos_plot <- function(data, contig_lengths,
                             tracks = c("Fst", "Snps_females", "Snps_males"),
                             track.labels = c("Fst", "SNP F.", "SNP M."),
                             track.colors = c("grey20", "firebrick2", "dodgerblue3"),
                             point.size = 0.1, bg.color = "white",
                             highlight = NULL, bg.highlight.color = "grey80",
                             output.file = NULL, width = 2400, height = 2400, res = 120,
                             sector.titles.expand = NULL) {

    # Check that sectors to highlight exist
    highlight <- check_highlight_sectors(data, contig_lengths, highlight)

    # Setup sector colors
    bg_colors <- rep(bg.color, length(contig_lengths))
    bg_colors[which(names(contig_lengths) %in% highlight)] <- bg.highlight.color

    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, length(contig_lengths)), contig_lengths), ncol=2)

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
    circlize::circos.initialize(factors = names(contig_lengths), xlim = sector_lengths, sector.width = contig_lengths)

    # Draw specified tracks
    top_track <- TRUE
    for (i in c(1:length(tracks))) {

        if (tracks[i] %in% names(data)) {

            track_data <- create_track_data(data, tracks[i], track.colors[i])
            track_circos_window(track_data, track.labels[i],
                                bg.col = bg_colors, point.size = 0.01,
                                top.track = top_track, sector.titles.expand = sector.titles.expand,
                                first_sector = names(contig_lengths)[1])

        } else {

            print(paste0("Warning: unknown track type \"", tracks[i], "\" ..."))

        }

        top_track <- FALSE
    }

    # Close output file
    if (!is.null(output.file)) { dev.off() }
}



check_highlight_sectors <- function(data, contig_lengths, highlight) {

    if (is.null(highlight)) {return(c())}

    for (i in 1:length(highlight)) {

        if (!(highlight[i] %in% names(contig_lengths))) {  # Sector to highlight not found in list of sector lengths

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



create_track_data <- function(data, track.name, track.color) {

    track_data <- data[, c("Contig_plot", "Position_plot", "Color", track.name)]
    track_data$Color = track.color

    return(track_data)
}