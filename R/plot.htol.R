#' @importFrom purrr map list_rbind
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @export
fortify.htol <- function(htol) {
  data <- map(htol, "data") |>
    list_rbind()

  preds <- map(htol, "preds") |>
    list_rbind(names_to = "id")

  crits <- map(htol, \(x) {
    tibble(
      T95 = mean(x$T95, na.rm = TRUE),
      T50 = mean(x$T50, na.rm = TRUE),
      Tcrit = mean(x$Tcrit, na.rm = TRUE)
    )
  }) |>
    list_rbind(names_to = "id") |>
    pivot_longer(-.data$id, names_to = "levels", values_to = "values")

  list(
    data = data,
    preds = preds,
    crits = crits
  )
}

#' @importFrom ggplot2 ggplot geom_point geom_ribbon aes facet_wrap vars
#' @export
autoplot.htol <- function(htol) {
  htol <- fortify.htol(htol)

  ggplot(htol$data, aes(x = .data$temperature)) +
    geom_point(aes(y = .data$fvfm)) +
    geom_ribbon(aes(ymax = .data$high, ymin = .data$low),
      data = htol$preds,
      alpha = 0.3, fill = "red"
    ) +
    geom_line(aes(y = .data$preds),
      data = htol$preds,
      colour = "red"
    ) +
    geom_vline(
      aes(xintercept = .data$values, colour = .data$levels),
      data = htol$crits,
      linetype = "dashed"
    ) +
    facet_wrap(vars(.data$id))
}
