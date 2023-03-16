#' A function to randomly select n colours and plot them.
#' Useful if you're stuck thinking of a palette
#' @title randomColours
#' @param number number of colours to pick
#' @param height number of rows for plot
#' @return a ggplot plot
#' @keywords colour
#' @export
#' @examples
#' randomColours()
#'
#' randomColours(number = 57, height = 8)


randomColours <- function(number = 10, height = 5) {
    colours <- data.frame("colourNumber" = 1:number,
                   colourCol = c(0:(number-1))%/%height,
                   colourRow = c(0:(number-1))%%height,
                   Colour = sample(colours(), number, replace = FALSE))
    ggplot(colours, aes(x = colourCol,
                        y = colourRow,
                        fill = Colour,
                        text = paste("<b>Colour:", Colour, "</b>"))) +
        geom_raster() +
        scale_fill_manual(values = setNames(colours$Colour, colours$Colour)) +
        shadowtext::geom_shadowtext(aes(label = Colour), color = "black",  fontface = "bold", bg.color = "white") +
        theme(legend.position = "none",
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_blank(),
              panel.grid = element_blank())
}
