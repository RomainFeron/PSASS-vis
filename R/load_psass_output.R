#' @title Load chromosome names
#'
#' @description Loads chromosome names from a tabulated file
#'
#' @param input_file_path Path to the chromosome names file
#'
#' @return A vector with contig IDs as names and chromosome names as values, or NULL if there is no input file
#'
#' @examples
#' chromosomes <- load_chromosome_names("chromosomes_names.tsv")

load_chromosome_names <- function(input_file_path) {

    if (!is.null(input_file_path)) {

        raw_data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE))
        data <- gtools::mixedsort(setNames(raw_data$X2, raw_data$X1))

    } else {

        data <- NULL

    }

    return(data)
}



#' @title Detect chromosomes
#'
#' @description Automatically detect chromosomes in a psass output data frame
#'
#' @param data Psass output data frame
#'
#' @return A named vector with contig IDs as names and chromosome names as values
#'
#' @examples
#' chromosomes <- detect_chromosomes(data)

detect_chromosomes <- function(data) {

    # Order contigs by length in a data frame
    contig_lengths <- data.frame(unique(data[, c(1, 3)]))[order(unique(data.frame(data[,3])), decreasing = TRUE),]

    # Identify chromosomes based on name: start with LG, CH, or NC (case-unsensitive)
    detected_chr <- subset(contig_lengths, toupper(substr(contig_lengths$Contig, 1, 2)) %in% c("LG", "CH", "NC"))

    if (nrow(detected_chr) > 0) {

        # Sometimes mitochondria are also called NC_xxx. If one chromosome is > 50 times smaller than the average of all other chromosomes,
        # or if it is smaller than 50000 bp, it is considered to be the mitochondrion and is removed
        potential_mt <- tail(detected_chr, 1)
        if (50 * potential_mt$Length < median(detected_chr$Length) | potential_mt$Length < 50000) detected_chr <- subset(detected_chr, detected_chr$Contig != potential_mt)

    }

    # Create the named vector of chromosome names with detected chromosomes
    return(setNames(detected_chr$Contig, detected_chr$Contig))
}



#' @title Load a genome input file
#'
#' @description Loads a file containing metrics along a genome. Format: Contig | Position | Length | <Metrics> ...
#'
#' @param input_file_path Path to a genome input file
#'
#' @param chromosomes Vector of chromosome names from \code{\link{load_chromosome_names}} (default: NULL)
#'
#' @param detect.chromosomes If TRUE, will consider contigs starting with LG, CH, or NC as chromosomes if no chromosomes were specified (default: TRUE)
#'
#' @return A list with two elements: "data" = parsed genome metrics data, "lengths" = lengths of contigs to be plotted
#'
#' @examples
#' window_data <- load_genome_input("psass_window.tsv", chromosomes=chromosomes)

load_genome_input <- function(input_file_path, chromosomes=NULL, detect.chromosomes=TRUE) {

    data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, trim_ws = TRUE))
    data$Contig_plot <- data$Contig
    data$Position_plot <- data$Position

    if (is.null(chromosomes) & detect.chromosomes) {

        chromosomes = detect_chromosomes(data)

    }

    if (!is.null(chromosomes)) {

        # Separate data between LGs and unplaced contigs
        lg <- subset(data, data$Contig %in% names(chromosomes))
        lg$Contig_plot <- chromosomes[lg$Contig_plot]
        unplaced <- subset(data, !(data$Contig %in% names(chromosomes)))

        # Set LG color index to 2 for plotting
        lg$Color <- rep(2, dim(lg)[1])

        # Initialize contig lengths vector with chromosome lengths
        plot_contig_lengths <- data.frame(unique(lg[, c(1, 3)]))
        plot_contig_lengths <- plot_contig_lengths[order(gtools::mixedorder(plot_contig_lengths$Contig)),]
        plot_contig_lengths$Contig <- chromosomes[plot_contig_lengths$Contig]

    } else {

        # No chromosomes: entire data is unplaced
        unplaced <- data
        lg <- NULL

    }

    # If there are unplaced contigs, concatenate these contigs in a chimeric scaffold
    if (nrow(unplaced) > 0) {

        # Order unplaced contigs data by contig length and then by position on the contig
        unplaced_len <- data.frame(unique(unplaced[, c(1, 3)]))[order(unique(data.frame(unplaced[,3])), decreasing = TRUE),]
        unplaced <- unplaced[order(match(unplaced$Contig, unplaced_len$Contig), unplaced$Position), ]

        # Attribute a color index to each unplaced contig, alternating between 0 and 1
        unplaced$Color <- setNames(seq(1, nrow(unplaced_len)), unplaced_len$Contig)[unplaced$Contig] %% 2

        # Transform position on each contig into position on cumulated contig
        cumulated_lengths <- setNames(c(0, head(cumsum(unplaced_len$Length), -1)), unplaced_len$Contig)
        unplaced$Position_plot <- unplaced$Position + cumulated_lengths[unplaced$Contig]
        unplaced$Contig_plot <- "Unplaced"

        # Compute total length of unplaced contigs
        total_unplaced_length <- sum(unplaced_len$Length)

    }

    # Create final data frame
    if (!is.null(lg) & nrow(unplaced) > 0) {

        data <- rbind(lg, unplaced)  # Combined lg and unplaced data
        plot_contig_lengths <- rbind(plot_contig_lengths, data.frame(Contig="Unplaced", Length=total_unplaced_length))

    } else if (!is.null(lg)) {

        data <- lg  # Entire dataset is in chromosomes

    } else if (nrow(unplaced) > 0) {

        data <- unplaced  # Entire dataset is unplaced
        plot_contig_lengths <- data.frame(Contig=c("Unplaced"), Length=c(total_unplaced_length))

    } else {

        print(paste0('Error loading dataset <', input_file_path, ">"))
        exit(1)

    }

    return(list(data=data, lengths=setNames(plot_contig_lengths$Length, plot_contig_lengths$Contig)))
}
