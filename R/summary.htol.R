#' Summary of htol objects
#' @description
#' Makes summary of htol objects made by `psiiht()`
#' @param object htol object make by `psiiht()`
#' @param \dots Extra arguments, currently unused

#' @export
summary.htol <- function(object, ...) {
  object$results |>
    map(\(x) {
      T50_ci <- quantile(x$T50, c(0.025, 0.975), na.rm = TRUE)


      data.frame(
        prop_success = mean(!is.na(x$T50)),
        Tcrit.mn = mean(na.omit(x$Tcrit)),
        T50.mn = mean(na.omit(x$T50)),
        T95.mn = mean(na.omit(x$T95))
      )
    }) |>
    list_rbind(names_to = "id")
}
