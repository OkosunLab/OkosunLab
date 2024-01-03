#' @title FancyUmapPlot
#' @param seuratObj a seurat object
#' @param title a title for the plot, if NULL a generic title will be added (default: NULL)
#' @param subtitle a subtitle for the plot, if NULL a generic subtitle will be added (default: NULL)
#' @param ... Any other arguments that will be passed onto dimplot
#' @return a ggplot object
#' @keywords UMAP dimplot
#' @export
#' @examples
#'
#' FancyUmapPlot(seuratObj)
#'
#' FancyUmapPlot(seuratObj, title = "test plot", subtitle = "a plot to test this function")
#'

FancyUmapPlot <-
    function (seuratObj,
              title = NULL,
              subtitle = NULL,
              ...) {
        umap <- Seurat::DimPlot(seuratObj, ...)
        if (! is.null(title)) {
            umap <- umap + labs(title = title)
        }
        if (! is.null(subtitle)) {
            umap <- umap + labs(subtitle = subtitle)
        }
        umap + theme(
            axis.line = element_line(arrow = arrow(
                length = unit(5,
                              "pt"), type = "closed"
            )),
            axis.text = element_blank(),
            axis.title = element_text(hjust = 0),
            axis.ticks = element_blank()
        ) +
            guides(
                x = ggh4x::guide_axis_truncated(
                    trunc_lower = unit(0,
                                       "npc"),
                    trunc_upper = unit(60, "pt")
                ),
                y = ggh4x::guide_axis_truncated(
                    trunc_lower = unit(0,
                                       "npc"),
                    trunc_upper = unit(60, "pt")
                )
            ) + labs(caption = paste("n =",
                                     nrow(seuratObj@meta.data), "cells"))
    }

#' @title FancyUmap
#' @param seuratObj a seurat object
#' @param split.by if Null proceed as normal (default). If set it will standardise the colours across the plots.
#' @param ... Any other arguments that will be passed onto dimplot
#' @return a ggplot object
#' @keywords UMAP dimplot
#' @export
#' @examples
#'
#' FancyUmapPlot(seuratObj)
#'
#' FancyUmapPlot(seuratObj, title = "test plot", subtitle = "a plot to test this function")
#'

FancyUmap <- function(seuratObj,
                      split.by = NULL,
                      ...) {
    if (is.null(split.by)) {
        umap <- FancyUmapPlot(seuratObj, ...)
    } else {
        args <- list(...)
        if (! "reduction" %in% names(args)) {
            reduction = "umap"
        } else {
            reduction = args$reduction
        }
        xlim <- Embeddings(seuratObj, reduction = reduction)[,1]
        ylim <- Embeddings(seuratObj, reduction = reduction)[,2]
        if ("cols" %in% names(args)) {
            umap <- lapply(unlist(unique(seuratObj[[split.by]])), function(split) {
                cells = colnames(seuratObj[,seuratObj[[split.by]] == split])
                FancyUmapPlot(seuratObj[,cells], ...) +
                    lims(x = c(floor(min(xlim)), ceiling(max(xlim))),
                         y = c(floor(min(ylim)), ceiling(max(ylim))))
            }) }
        else {
            if (! "group.by" %in% names(args)){
                groups = unlist(unique(seuratObj@active.ident))
            } else {
                groups = unlist(unique(seuratObj[[args$group.by]]))
            }
            ## If you don't specific the colours you will add NAs here which will get a colour
            ## Remove them here
            groups <- groups[!is.na(groups)]
            if (! "na.value" %in% names(args)) {
                na.value = "gray50"
            } else {
                na.value = args$na.value
            }
            colours <- scales::hue_pal()(length(groups)) %>%
                setNames(groups)
            umap <- lapply(unlist(unique(seuratObj[[split.by]])), function(split) {
                cells = colnames(seuratObj[,seuratObj[[split.by]] == split])
                FancyUmapPlot(seuratObj[,cells], ...) +
                    lims(x = c(floor(min(xlim)), ceiling(max(xlim))),
                         y = c(floor(min(ylim)), ceiling(max(ylim)))) +
                    scale_colour_manual(values = colours,
                                        drop = FALSE,
                                        na.value = na.value)
            })}
        umap <- umap %>% patchwork::wrap_plots()
    }
    umap
}
