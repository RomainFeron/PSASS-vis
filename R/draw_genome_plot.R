
check_highlight_sectors <- function(data, contig_lengths, highlight) {

    # Test if lgs specified in "highlight" exist, otherwise remove them from the highlight vector and print a warning
    to_remove <- c()

    for (i in 1:length(highlight)) {

        if (!(highlight[i] %in% names(contig_lengths))) {  # Sector to highlight not found in list of sector lengths

            contigs <- setNames(unique(data$Contig_plot), unique(data$Contig))

            if (highlight[i] %in% names(contigs)) {  # Look for specified sector to highlight in original contig names

                highlight[i] <- contigs[highlight[i]]

            } else {

                warning(paste0("Could not find sector to highlight \"", highlight[i], "\"."))
                to_remove <- c(to_remove, i)

            }
        }
    }

    if (!is.null(to_remove)) {highlight <- highlight[-to_remove]}  # Renmove sectors not found in new and old contig names

    return(highlight)
}


update_data_with_highlights <- function(data, highlight, suffix=NULL) {

    highlight_data <- data[data$Contig_plot %in% highlight, , drop = FALSE]  # Extract highlighted sectors data
    highlight_data$Contig_plot <- paste0(highlight_data$Contig_plot, suffix)  # Add suffix to highlighted contig names
    data <- rbind(data, highlight_data)  # Combine base sectors data and highlighted sectors data

    return(data)
}


update_lengths_with_highlights <- function(contig_lengths, highlight, suffix=NULL) {

    highlight_lengths <- contig_lengths[highlight]  # Extract highlighted sector lengths
    names(highlight_lengths) <- paste0(highlight, suffix)  # Edit highlighted sector names with suffix
    contig_lengths <- c(contig_lengths, highlight_lengths)  # Combine base sectors lengths and highlighted sectors lengths

    return(contig_lengths)
}



data = psass_window_chr$data
contig_lengths = psass_window_chr$lengths
output.file = NULL
width = 2400
height = 2400
res = 120
highlight = c("Chr01", "Chr10")
expand.highlight.suffix = " (highlight)"
expand.highlight.ratio = 1
bg.color = "white"
bg.highlight.color = "grey80"
point.size = 0.1
points.color = "grey20"
pools.color = c("firebrick2", "dodgerblue3")
color.unplaced = FALSE
points.color.unplaced = c("dodgerblue3", "goldenrod1")
sector.titles.expand = NULL
tracks = c("fst", "snp_males")


draw_circos_plot <- function(data, contig_lengths,
                             output.file = NULL, width = 2400, height = 2400, res = 120,
                             tracks = c("fst_win", "snp_win_pool1", "snp_win_pool2", "depth_ratio"),
                             highlight = NULL, expand.highlight.ratio = 0, expand.highlight.suffix = " (highlight)",
                             bg.color = "white", bg.highlight.color = "grey80", point.size = 0.1,
                             points.color = "grey20",
                             pools.color = c("firebrick2", "dodgerblue3"),
                             color.unplaced = FALSE, points.color.unplaced = c("dodgerblue3", "goldenrod1"),
                             sector.titles.expand = NULL) {

    # Check that sectors to highlight exist
    highlight <- check_highlight_sectors(data, contig_lengths, highlight)

    if (expand.highlight.ratio > 0) {
        # Update data with highlighted sectors
        data <- update_data_with_highlights(data, highlight, suffix=expand.highlight.suffix)
        # Update data with highlighted sectors
        contig_lengths <- update_lengths_with_highlights(contig_lengths, highlight, suffix=expand.highlight.suffix)
    }

    # Setup sector colors
    bg_colors <- rep(bg.color, length(contig_lengths))
    bg_colors[which(names(contig_lengths) %in% highlight)] <- bg.highlight.color
    bg_colors[which(names(contig_lengths) %in% paste0(highlight, expand.highlight.suffix))] <- bg.highlight.color

    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, length(contig_lengths)), contig_lengths), ncol=2)

    # Setup gaps between sectors
    gaps <- rep(1, nrow(sector_lengths))
    gaps[length(gaps)] <- 12  # Bigger gap after the last sector to leave some space for track titles
    if (expand.highlight.ratio > 0 & (length(highlight) > 0)) { gaps[length(gaps) - length(highlight)] <- 5 }  # Add a small gap before the first zoomed sector

    # Calculate angle offset to have zoomed sector on top
    if (expand.highlight.ratio > 0 & length(highlight) > 0) {
        a <- 360 * sum(tail(contig_lengths, n = length(highlight))) / sum(contig_lengths) / 2 + 7.5
    } else {
        a <- 7.5
    }

    # Calculate width of track and sector.title.expand factor based on number of tracks
    n_tracks <- length(tracks)
    track.height <- 0.25  # Width of 0.25 if single track is plotted
    if (n_tracks >= 2) { track.height <- 0.5 / n_tracks }  # For multiple tracks, total width 0.5 (half the circo's width)
    if (is.null(sector.titles.expand)) sector.titles.expand <- 1.1 + n_tracks / 10

    # Reset circos parameters
    circlize::circos.clear()

    # Setup circos parameters
    circlize::circos.par("track.height" = track.height,
                         "start.degree" = 90 - a,
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
    top.track <- TRUE
    for (i in c(1:length(tracks))) {

        if (TRUE) {


        } else {

            print(paste0("Warning: unknown track type \"", tracks[i], "\" ..."))

        }

        top.track <- FALSE
    }

    # Close output file
    if (!is.null(output.file)) {
        dev.off()
    }
}
