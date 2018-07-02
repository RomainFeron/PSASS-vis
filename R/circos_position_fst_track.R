#' @title Circos FST position track
#'
#' @description Draws a track on a circos plot for FST position data. This function is intended for use in the \code{\link{circos_plot}} function.
#'
#' @param data FST position data frame.
#'
#' @param ylim Limits of the y axis (default c(min(FST), 1)).
#'
#' @param bg.col Background color for sectors, either a single color or a vector of colors for each sector (default "white").
#'
#' @param point.size Size of a point in the plot (default 0.5).
#'
#' @param top.track If TRUE, this track will be considered top track of the plot and the x axis will be drawn (default FALSE).
#'
#' @param sector.names Vector of contig names obtained with the \code{\link{load_contig_names}} (default NULL).
#'
#' @param sector.title.expand Value controlling the distance between sector titles and the top axis (default 1.3).
#'
#' @param color.palette Color palette for this track (default c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey20")).


draw_position_fst <- function(data,
                              ylim = c(min(data$Fst), 1), bg.col = "white", point.size = 0.5,
                              top.track = FALSE, sector.names = NULL, sector.titles.expand = 1.3,
                              color.palette = c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey20") ) {

    # Draw the top track of the plot, showing sex-bias
    circlize::circos.track(factors = data$Contig,
                           x = data$Position,
                           y = data$Fst,
                           ylim = ylim,
                           bg.col = bg.col,
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
                                   if (!is.null(sector.names)) {
                                       circlize::circos.text(xcenter,
                                                             sector.titles.expand * ymax,
                                                             sector.names[sector.index],
                                                             cex = 1.5,
                                                             facing = "bending.inside",
                                                             niceFacing = TRUE)
                                   }
                               }

                               # Plot the data
                               circlize::circos.points(x, y, cex = point.size,
                                                       col = color.palette[as.character(data$Color[which(data$Contig == sector.index)])],
                                                       bg = color.palette[as.character(data$Color[which(data$Contig == sector.index)])],
                                                       pch = 21)

                               # Add Y axis on the first sector only
                               if (sector.index == sectors[1]) {

                                   # Create y axis
                                   circlize::circos.yaxis(side = "left",
                                                          at = c(ylim[1], (ylim[2] - ylim[1]) / 2 + ylim[1], ylim[2]),  # 3 labels
                                                          sector.index = sectors[1],
                                                          labels.cex = 1.2,
                                                          labels.niceFacing = FALSE,
                                                          labels = round(c(ylim[1], (ylim[2] - ylim[1]) / 2 + ylim[1], ylim[2]), 2))

                                   #Add y axis labels
                                   label_offset <- - 5 * (xmax - xmin) / (xplot[1] - xplot[2])  # Axis title will be plotted 5° on the left of the axis
                                   circlize::circos.text(label_offset,
                                                         0.5 * (ymax - ymin) + ymin,
                                                         expression(bold(paste0("F"["ST"], " > 0.25"))),
                                                         sector.index = sectors[1],
                                                         facing = "clockwise",
                                                         cex = 1.3,
                                                         font = 2)
                               }
                           }
    )

}