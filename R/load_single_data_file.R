#' @title Load psass sliding window output
#'
#' @description Loads the main output file from psass analyze function (sliding window metrics)
#'
#' @param input_file_path Path to the output of psass
#'
#' @param chromosomes Vector of chromosome names from \code{\link{load_chromosome_names}}
#'
#' @return A data frame storing the data from the input file
#'
#' @examples
#' window_data <- load_psass_window_output("psass_window.tsv")


input_file_path = "psass_window.tsv"

load_psass_window_output <- function(input_file_path, chromosomes=NULL) {

    data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, trim_ws = TRUE))

    if (chromosomes) {

        # Separate data between LGs and unplaced contigs
        lg <- subset(data, data$Contig %in% names(chromosomes))
        unplaced <- subset(data, !(data$Contig %in% names(chromosomes)))

        # Set LG color index to 2 for plotting
        lg$Color_index <- rep(2, dim(lg)[1])

        if (nrow(unplaced) > 0) {

            # Order unplaced contigs data by contig length and then by position on the contig
            unplaced <- unplaced[order(match(unplaced$Contig, names(contig_lengths$unplaced)), unplaced$Position), ]
        }

        # Attribute a color index to each unplaced contig, alternating between 0 and 1
        order <- seq(1, length(unique(data_unplaced$Contig)))
        names(order) <- unique(data_unplaced$Contig)
        data_unplaced$Color <- order[data_unplaced$Contig] %% 2

        # Transform position on each contig into position on cumulated contig
        temp <- cumsum(contig_lengths$unplaced) - contig_lengths$unplaced[1]
        data_unplaced$Original_position <- data_unplaced$Position
        data_unplaced$Position <- data_unplaced$Position + temp[data_unplaced$Contig]
        data_unplaced$Contig_id <- data_unplaced$Contig
        data_unplaced$Contig <- "Unplaced"

        # Regroup data into one data frame
        lg$Contig_id <- lg$Contig
        lg$Original_position <- lg$Position
        data <- rbind(lg, data_unplaced)
        data$Contig <- factor(data$Contig, levels = c(names(contig_lengths$lg), "Unplaced"))
    }

    if (plot.unplaced & dim(data_unplaced)[1] > 0) {  # If unplaced scaffolds should be grouped and there is at least one unplaced scaffold




        # Attribute a color index to each unplaced contig, alternating between 0 and 1
        order <- seq(1, length(unique(data_unplaced$Contig)))
        names(order) <- unique(data_unplaced$Contig)
        data_unplaced$Color <- order[data_unplaced$Contig] %% 2

        # Transform position on each contig into position on cumulated contig
        temp <- cumsum(contig_lengths$unplaced) - contig_lengths$unplaced[1]
        data_unplaced$Original_position <- data_unplaced$Position
        data_unplaced$Position <- data_unplaced$Position + temp[data_unplaced$Contig]
        data_unplaced$Contig_id <- data_unplaced$Contig
        data_unplaced$Contig <- "Unplaced"

        # Regroup data into one data frame
        lg$Contig_id <- lg$Contig
        lg$Original_position <- lg$Position
        data <- rbind(lg, data_unplaced)
        data$Contig <- factor(data$Contig, levels = c(names(contig_lengths$lg), "Unplaced"))

    } else {

        data <- lg
        data$Original_position <- data$Position
        data$Contig <- factor(data$Contig, levels = names(contig_lengths$lg))

    }

    return(data)
}
