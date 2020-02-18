#' @title Plot a circos track
#'
#' @description Plot a single track for a circos plot
#'
#' @param data Input data generated with the \code{\link{create_region_track_data}} function)
#'
#' @param track Track object storing properties for the current track, generated with the \code{\link{circos_track}} function
#'
#' @param top.track If TRUE, x-axis labels and sector names will be added to the track
#'
#' @param sector.titles.expand Manually set the space between sector titles and x-axis as a multiple of ymax (default: 1.3)
#'
#' @param first.sector Name of the first sector to draw y-axis on (default: NULL)
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' fst_track <- circos_track("Fst")
#' circos_data <- create_circos_track_data(genomic_data, fst_track)
#'
#' plot_track_region(region_data, fst_track, top.track=TRUE)
#'

plot_circos_track <- function(data, track, top.track = FALSE, sector.titles.expand = 1.3, first.sector=NULL) {

    # Assign values for y-axis parameters
    if (is.null(track$ylim)) { track$ylim = c(0.975 * min(data[, 4]) - 0.01, 1.025 * max(data[, 4]) + 0.01) }
    ylim <- track$ylim
    ylabel <- track$label

    # Draw the top track of the plot
    circlize::circos.track(factors = data$Contig_plot,
                           x = data$Position_plot,
                           y = as.vector(unlist(data[,4])),
                           ylim = track$ylim,
                           bg.col = track$bg_colors,
                           panel.fun = function(x, y) {  # panel.fun is the function drawing the track

                               # Get useful sector information
                               sector.index <- circlize::get.cell.meta.data("sector.index")
                               xcenter <- circlize::get.cell.meta.data("xcenter")
                               ymin <- circlize::get.cell.meta.data("ylim")[1]
                               ymax <- circlize::get.cell.meta.data("ylim")[2]
                               xmin <- circlize::get.cell.meta.data("xlim")[1]
                               xmax <- circlize::get.cell.meta.data("xlim")[2]
                               xplot <- circlize::get.cell.meta.data("xplot")

                               # Add top axis and titles to sectors
                               if (top.track) {

                                   # Create x axis on top of sectors
                                   circlize::circos.axis(h = "top",
                                                         major.at = c(0, xmax / 3, 2 * xmax / 3, xmax),  # Label every 1/3 of the axis
                                                         labels.cex = 1.2,
                                                         labels.facing = "outside",
                                                         direction="outside",
                                                         labels = convert_to_mb(c(0, xmax / 3, 2 * xmax / 3, xmax)),  # Conversion to Mb
                                                         minor.ticks = 4,
                                                         labels.pos.adjust = TRUE)

                                   # Add sector names
                                   circlize::circos.text(xcenter,
                                                         sector.titles.expand * ymax,
                                                         sector.index,
                                                         cex = 1.5,
                                                         facing = "bending.inside",
                                                         niceFacing = TRUE)
                               }

                               # Plot the data
                               circlize::circos.points(x, y, cex = track$point.size,
                                                       col = data$Color,
                                                       bg = data$Color,
                                                       pch = 21)

                               # Add Y axis on the first sector only
                               if (sector.index == first.sector) {

                                   # Create y axis
                                   circlize::circos.yaxis(side = "left",
                                                          at = c(ylim[1], (ylim[2] - ylim[1]) / 2 + ylim[1], ylim[2]),  # 3 labels
                                                          sector.index = first.sector,
                                                          labels.cex = 1.2,
                                                          labels.niceFacing = FALSE,
                                                          labels = round(c(ylim[1], (ylim[2] - ylim[1]) / 2 + ylim[1], ylim[2]), 0))

                                   #Add y axis labels
                                   label_offset <- - 7.5 * (xmax - xmin) / (xplot[1] - xplot[2])  # Axis title will be plotted 7.5° on the left of the axis
                                   circlize::circos.text(label_offset,
                                                         0.5 * (ymax - ymin) + ymin,
                                                         ylabel,
                                                         sector.index = first.sector,
                                                         facing = "clockwise",
                                                         cex = 1.3,
                                                         font = 2)
                               }
                           }
    )
}