#' @title Load a chromosomes names file
#'
#' @description Loads a table of chromosomes names.
#'
#' @param input_file_path Path to a chromosomes names file.
#'
#' @return A vector with contigs as names and chromosomes names as values, or NULL if there is no input file
#'
#' @examples
#' lengths <- load_chromosomes_names("chromosomes_names.tsv")



load_chromosomes_names <- function(input_file_path) {

    if (!is.null(input_file_path)) {

        raw_data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE))
        data <- raw_data$X2
        names(data) <- raw_data$X1
        data <- gtools::mixedsort(data)
        return(data)

    } else {

        return(NULL)

    }
}

