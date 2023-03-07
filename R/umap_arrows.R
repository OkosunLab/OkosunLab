#' @title umap_arrows
#' @param plot a ggplot object
#' @param percentage a value between 0 and 1 for how far on each axis the arrows should reach
#' @param label_x a string for the x arrow text
#' @param label_y a string for the y arrow text
#' @param x_pad padding for the x arrow label (percentage of arrow length)
#' @param y_pad padding for the y arrow label (percentage of arrow length)
#' @return a ggplot object
#' @keywords UMAP dimplot
#' @export
#' @examples
#' p <- ggplot(mtcars, aes(x = mpg, y = wt)) +
#'     geom_point()
#' umap_arrows(p)

umap_arrows <- function(plot, percentage = 0.2, label_x = "UMAP 1", label_y = "UMAP 2", x_pad = 0.5, y_pad = 0.5) {
    xrange <- layer_scales(plot)$x$range$range
    yrange <- layer_scales(plot)$y$range$range

    x_wid <-
        percentage * (xrange[2] - xrange[1])
    y_wid <-
        percentage * (yrange[2] - yrange[1])

    x_pad <- x_wid * x_pad
    y_pad <- y_wid * y_pad

    xstart <- layer_scales(plot)$x$range$range[1]
    xstop <- xstart + x_wid
    ystart <- layer_scales(plot)$y$range$range[1]
    ystop <- ystart + y_wid

    plot +
        annotate(
            "segment",
            x = xstart,
            xend = xstop + c(0,-x_wid),
            y = ystart,
            yend = ystop + c(-y_wid, 0),
            arrow = arrow(type = "closed", length = unit(10, 'pt'))
        ) +
        annotate(
            "text",
            x = xstart + x_wid / 2,
            y = ystart - y_pad,
            label = "label_x",
            size = 3,
            vjust = 0.5,
            hjust = 0.5
        ) +
        annotate(
            "text",
            x = xstart - x_pad,
            y = ystart + y_wid / 2,
            label = "label_y",
            angle = 90,
            size = 3,
            vjust = 0.5,
            hjust = 0.5
        )
}
