#' Linköping University ggplot2 theme
#' @return A ggplot2 theme object
#' @export
liu_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#d1f4f6"),

      plot.title = ggplot2::element_text(face = "bold"),

      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#54d8e0"),
    )
}
