#' @title umap_arrows
#' @param seuratObj a seurat object
#' @param reduction the reduction to use for the Dimplot (default: umap)
#' @param group the column in the metadata to colour the cells by (default: orig.ident)
#' @param title a title for the plot, if NULL a generic title will be added (default: NULL)
#' @param subtitle a subtitle for the plot, if NULL a generic subtitle will be added (default: NULL)
#' @param split a column in the metadata to split the umap into multiple umaps (default: NULL) N.B. this will break the arrows!
#' @return a ggplot object
#' @keywords UMAP dimplot
#' @export
#' @examples
#'
#' FancyUmap(seuratObj, reduction = "umap", group = "seurat_clusters", title = "test plot", subtitle = "a plot to test this function", split = NULL)
#'


FancyUmap <- function(seuratObj,
                      reduction = "umap",
                      group = "orig.ident",
                      title = NULL,
                      subtitle = NULL,
                      split = NULL) {
    if (is.null(split)) {
        umap <- DimPlot(seuratObj,
                        reduction = "umap", group.by = group)
    } else {
        umap <- DimPlot(seuratObj,
                        reduction = "umap",
                        group.by = group,
                        split.by = split)
    }
    if (is.null(title)) {
        umap <- umap +
            labs(title = paste("A", reduction, "plot"))
    } else {
        umap <- umap +
            labs(title = title)
    }
    if (is.null(subtitle)) {
        umap <- umap +
            labs(subtitle = paste("Cells coloured by", group))
    } else {
        umap <- umap +
            labs(subtitle = subtitle)
    }
    ## Completely remove axes
    umap +
        theme(
            axis.line = element_line(arrow = arrow(
                length = unit(5, "pt"),
                type = "closed"
            )),
            axis.text = element_blank(),
            axis.title = element_text(hjust = 0),
            axis.ticks = element_blank()
        ) +
        ## Use ggh4x's truncated axis to plot only 60pt axes
        guides(
            x = ggh4x::guide_axis_truncated(
                trunc_lower = unit(0, "npc"),
                trunc_upper = unit(60, "pt")
            ),
            y = ggh4x::guide_axis_truncated(
                trunc_lower = unit(0, "npc"),
                trunc_upper = unit(60, "pt")
            )
        ) +
        labs(caption = paste("n =",
                             nrow(seuratObj@meta.data),
                             "cells"))
}
