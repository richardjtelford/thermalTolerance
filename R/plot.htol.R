#' Fortify htol objects
#' @description
#' Process htol objects made by `psiiht()` for `autoplot()`
#' @param model htol object make by `psiiht()`
#' @param data Orginal data - provided by default
#' @param \dots Extra arguments, currently unused

#' @importFrom purrr map list_rbind
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 fortify
#' @importFrom rlang .data
#' @export
fortify.htol <- function(model, data = model$data, ...) {


  preds <- map(model$results, "preds") |>
    list_rbind(names_to = "id")

  crits <- map(model$results, \(x) {
    tibble(
      T95 = mean(x$T95, na.rm = TRUE),
      T50 = mean(x$T50, na.rm = TRUE),
      Tcrit = mean(x$Tcrit, na.rm = TRUE)
    )
  }) |>
    list_rbind(names_to = "id") |>
    pivot_longer(-.data$id, names_to = "levels", values_to = "values")

  list(
    data = model$data,
    preds = preds,
    crits = crits
  )
}

#' Autoplot htol objects
#' @description
#' Default plots for htol objects made by `psiiht`
#' 
#' @param object Result from `psiiht`
#' @param \dots Extra arguments, currently unused 
#' @importFrom ggplot2 ggplot geom_point geom_ribbon aes facet_wrap vars geom_line geom_vline autoplot
#' @export
autoplot.htol <- function(object, ...) {
  object <- fortify.htol(object)

  ggplot(object$data, aes(x = .data$temperature)) +
    geom_point(aes(y = .data$fvfm)) +
    geom_ribbon(aes(ymax = .data$high, ymin = .data$low),
      data = object$preds,
      alpha = 0.3, fill = "red"
    ) +
    geom_line(aes(y = .data$preds),
      data = object$preds,
      colour = "red"
    ) +
    geom_vline(
      aes(xintercept = .data$values, colour = .data$levels),
      data = object$crits,
      linetype = "dashed"
    ) +
    facet_wrap(vars(.data$id))
}
