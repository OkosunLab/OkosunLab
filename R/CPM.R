#' A function to generate cpm from expression matrix
#'
#'
#' @title CPM
#' @param matrix a matrix of counts, rows for features and columns for samples
#' @param log boolean whether to log2 transform the cpm
#' @param pseudocount boolean whether to add a pseudocount of 1 to the matrix
#' @return a normalised count matrix
#' @keywords normalisation cpm expression
#' @export
#' @examples
#' CPM(mtcars)
#'
#' CPM(mtcars, log = TRUE)
#'
#' CPM(mtcars, log = TRUE, pseudocount = TRUE)

CPM <- function(matrix, log = FALSE, pseudocount = FALSE) {
    ## You can't log with zeroes so if you have zeroes you can add a pseudocount of 1 here
    if ( pseudocount ) {
        matrix <- matrix + 1
    }
    ## Calculate the per million scaling factor
    sf <- colSums(matrix)/(1*10^6)
    ## divide each row by the vector of scaling factors
    ## returns rows as columns
    cpm <- apply(matrix, 1, function(col) { col/sf })
    ## transpose the cpm back to normal
    cpm <- t(cpm)
    ## If log is true run log on the CPM
    if ( log ) {
        cpm <- log2(cpm)
    }
    ## return cpm
    cpm
}
