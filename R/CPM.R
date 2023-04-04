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


#' A function to generate tpm from expression matrix
#'
#'
#' @title TPM
#' @param matrix a matrix of counts, rows for features and columns for samples
#' @param length a vector of gene lengths in the same order as the matrix
#' @param log boolean whether to log2 transform the cpm
#' @param pseudocount boolean whether to add a pseudocount of 1 to the matrix
#' @return a normalised count matrix
#' @keywords normalisation tpm expression
#' @export
#' @examples
#' TPM(matrix, lengths)
#'

TPM <- function(matrix, length, log = TRUE, pseudocount = TRUE) {
    ## You can't log with zeroes so if you have zeroes you can add a pseudocount of 1 here
    if ( pseudocount ) {
        matrix <- matrix + 1
    }
    ## reads / gene length
    rpk <- matrix/length
    ## scaling factors
    sf <- colSums(rpk)/(1*10^6)
    ## divide rpk by sf
    ## returns feature as column and sample as row so transpose
    tpm <- apply(rpk, 1, function(col) {col/sf})
    ## returns feature as column and sample as row so transpose
    tpm <- t(tpm)
    if (log) {
        log2(tpm)
    }
}
