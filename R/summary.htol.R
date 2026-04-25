#' Summary of htol objects
#' @description
#' Makes summary of htol objects made by `psiiht()`
#' @param object htol object make by `psiiht()`
#' @param \dots Extra arguments, currently unused

#' @export
summary.htol <- function(object, ...) {
  object$results |>
    map(\(x) {
      T50_ci <- quantile(x$T50, c(0.025, 0.5, 0.975), na.rm = TRUE) |> 
        unname()


      data.frame(
        prop_success = mean(!is.na(x$T50)),
        Tcrit.mn = mean(x$Tcrit, na.rm = TRUE),
        T50.mn = mean(x$T50, na.rm = TRUE),
        T50_ci_low = T50_ci[1],
        T50_ci_median = T50_ci[2],
        T50_ci_high = T50_ci[3],
        T95.mn = mean(x$T95, na.rm = TRUE)
      )
    }) |>
    list_rbind(names_to = "id")
}
