#' project values to a Desikan-Killiany brain atlas
#' 
#' @param dat A data frame.
#' @param fill A string.
#' @param title A string.
#' @param min A numeric.
#' @param max A numeric.
#' @returns A ggplot figure.
#' 
#' 
#' @importFrom magrittr "%>%"

library(ggseg)

plot_ggseg_brain <- function(
    dat, atlas = "dk", fill, min, max, break_int, title = NULL, 
    fill_label = NULL
) {
    # color bar set up
    colfunc <- colorRampPalette(c("#395D9C", "#358CA7", "white", "#F57A17", "#CE204E"))
    
    colors <- colfunc(100)

    if (atlas == "aseg") {
        fig <- dat %>%
            ggplot() +
            geom_brain(
                atlas = aseg, 
                position = position_brain("coronal"),
                mapping = aes(
                    fill = !!sym(fill)
                )
            ) +
            theme_void()
    } else if (atlas == "dk") {
        fig <- dat %>%
            ggplot() +
            geom_brain(
                atlas = dk, 
                position = position_brain(. ~ hemi + side),
                mapping = aes(
                    fill = !!sym(fill)
                )
            ) +
            scale_x_continuous(
                breaks = c(350, 1050),
                labels = c("Left", "Right")
            ) +
            theme_void() +
            theme(
                axis.text.x = element_text(),
                legend.key.size = unit(5, "mm")
            )
    }
    
    fig <- fig +
        scale_fill_gradientn(
            colors = colors,
            na.value = "grey85",
            limits = c(min, max),
            breaks = seq(min, max, break_int)
        ) +
        labs(title = title, fill = fill_label) +
        theme(
            plot.margin = margin(10, 10, 10, 10, "mm")
        )
    
    return(fig)
}
