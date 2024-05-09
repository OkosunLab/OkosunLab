#' A function to Calculate QC metrics from a list  of seurat objects
#' Useful if you're stuck thinking of a palette
#'
#' @title CalculateSeuratQC
#' @param seuratList a list of seurat objects
#' @return a list of seurat objects
#' @keywords scRNAseq
#' @export
#' @examples
#' CalculateSeuratQC()
#'
#' CalculateSeuratQC(seuList)
#'

CalculateSeuratQC <- function(seuratObj) {
    # Add number of genes per UMI for each cell to metadata
    seuratObj$log10GenesPerUMI <-
        log10(seuratObj$nFeature_RNA) / log10(seuratObj$nCount_RNA)
    ## Get percentage of genes starting MT- (mitochondrial genes in humans)
    seuratObj[["percent.mt"]] <-
        PercentageFeatureSet(seuratObj, pattern = "^MT-")
    ## Return to ratio
    seuratObj[["mitoRatio"]] <-
        seuratObj[["percent.mt"]] / 100
    ### Ribosomal ratio
    ## Get percentage of ribosomal genes
    seuratObj[["percent.ribo"]] <-
        PercentageFeatureSet(seuratObj, pattern = "^RP[S,L]")
    ## Return to ratio
    seuratObj[["riboRatio"]] <-
        seuratObj[["percent.ribo"]] / 100
    seuratObj
}


#' A function to plot QC metrics from a list  of seurat objects
#' Useful if you're stuck thinking of a palette
#'
#' @title PlotSeuratQCFromList
#' @param seuratList a list of seurat objects
#' @return a ggplot plot
#' @keywords scRNAseq
#' @export
#' @examples
#' PlotSeuratQCFromList()
#'
#' PlotSeuratQCFromList(seuList)

PlotSeuratQCFromList <- function(seuratList) {
    ## Get a metadata df for QC plots
    combined.meta.data <- lapply(seuratObjs,
                                 function(seuratObj) {
                                     seuratObj@meta.data
                                 }) %>%
        bind_rows() %>%
        mutate(log10nFeature_RNA = log10(nFeature_RNA),
               log10nCount_RNA = log10(nCount_RNA))
    ## These are the cut offs I plan to use below
    MetricCuttofs <-
        data.frame(
            Metrics = c("log10nCount_RNA", "log10nFeature_RNA", "log10GenesPerUMI", "mitoRatio"),
            Cutoffs = c(log10(200), log10(200), 0.8, 0.1)
        )
    combined.meta.data %>%
        pivot_longer(
            c(
                log10nCount_RNA,
                log10nFeature_RNA,
                log10GenesPerUMI,
                mitoRatio,
                riboRatio
            ),
            names_to = "Metrics",
            values_to = "Value"
        ) %>%
        ggplot(aes(x = Value,
                   y = interaction(Patient, Tissue, sep = "&"),
                   fill = Patient)) +
        ggridges::geom_density_ridges(alpha = 0.5) +
        facet_grid(Sort ~
                       fct_relevel(Metrics,
                                   c(
                                       "log10nCount_RNA",
                                       "log10nFeature_RNA",
                                       "log10GenesPerUMI",
                                       "mitoRatio",
                                       "riboRatio"
                                   )), scales = "free_x") +
        geom_vline(
            data = MetricCuttofs,
            aes(xintercept = Cutoffs),
            linetype = "dashed",
            colour = "red"
        ) +
        guides(y = ggh4x::guide_axis_nested(delim = "&")) +
        theme_classic() +
        labs(y = "Sample")
}
