#' Visualise Average Flight Delays by Airport Location
#'
#' Aggregates the `nycflights13::flights` data set using `dplyr` and returns a
#' `ggplot2` object showing the mean arrival delay by airport longitude and
#' latitude. Point size encodes the number of flights contributing to the
#' average. Adds a U.S. country outline in the background.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @importFrom dplyr filter group_by summarise n inner_join mutate arrange desc
#' @importFrom stats reorder
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
visualize_airport_delays <- function() {
  if (!requireNamespace("nycflights13", quietly = TRUE)) {
    stop("Package 'nycflights13' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("maps", quietly = TRUE)) {
    stop("Package 'maps' is required for the background outline. Install it with install.packages('maps').",
         call. = FALSE)
  }

  airports <- nycflights13::airports
  flights  <- nycflights13::flights

  delay_summary <- flights |>
    filter(!is.na(arr_delay)) |>
    group_by(dest) |>
    summarise(
      mean_delay   = mean(arr_delay, na.rm = TRUE),
      flight_count = n(),
      .groups = "drop"
    ) |>
    inner_join(airports, by = c("dest" = "faa")) |>
    mutate(dest = reorder(dest, mean_delay)) |>
    arrange(desc(mean_delay))

  ggplot2::ggplot(delay_summary, aes(x = .data$lon, y = .data$lat)) +
    ggplot2::borders("state", fill = NA, colour = "grey80") +
    ggplot2::geom_point(
      aes(colour = .data$mean_delay, size = .data$flight_count),
      alpha = 0.8
    ) +
    ggplot2::scale_colour_gradient2(
      low = "steelblue",
      mid = "khaki",
      high = "firebrick",
      midpoint = 0, # center at on-time
      name = "Mean arr. delay"
    ) +
    ggplot2::scale_size_continuous(name = "Flights", range = c(1, 6)) +
    ggplot2::labs(
      title = "Average Arrival Delay by Destination",
      subtitle = "NYC flights in 2013",
      x = "Longitude",
      y = "Latitude"
    ) +
    ggplot2::coord_quickmap(expand = FALSE) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1)) +
    ggplot2::theme_minimal()
}

visualize_airport_delays()
