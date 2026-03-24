#' @importFrom purrr map list_rbind
#' @export

fortify.htol <- function(htol) {
  print("in fortify")
  data <- map(htol, "data") |> 
    list_rbind()
  
  fvfm.boot <- map(htol, "fvfm.boot") |> 
    list_rbind(names_to = "id")
  
  list(
    data = data, 
    fvfm.boot = fvfm.boot
  )
}

#' @importFrom ggplot2 ggplot geom_point geom_ribbon aes facet_wrap vars
#' @export
autoplot.htol <- function(htol) {
  htol <- fortify.htol(htol)
  
  ggplot(htol$data, aes(x = temperature)) +
    geom_point(aes(y = fvfm)) +
    geom_ribbon(aes(ymax = high, ymin = low), data = htol$fvfm.boot, 
                alpha = 0.3, fill = "red") +
    facet_wrap(vars(id))
}

plot.htol <- function(htol) {

  
  plot(NULL, NULL, xlab = "Temperature", ylab = "Fv/Fm", xlim = c(23, 65), ylim = c(0, 0.9), bty = "l", lty = 2)
  text(22, 0.3, pos = 4, paste(unique(paste(id))), font = 4)
  points(Temperature, fvfm, xlab = "Temperature", ylab = "Fv/Fm", xlim = c(23, 65), ylim = c(0, 0.85), bty = "l", lty = 2, pch = 5, col = "black")
  lines(seq(23, 62), y, lwd = 1, col = "black")
  if (boots > 1) {
    lines(seq(23, 62), fvfm.boot[, 1], lty = 3, col = "black")
    lines(seq(23, 62), fvfm.boot[, 2], lty = 3, col = "black")
  }
  text(30, 0, paste("Tcrit:", round(mean(na.omit(Tcrit)), 1)), pos = 4, col = "light gray")
  abline(v = round(mean(na.omit(Tcrit)), 1), lty = 2, lwd = 1.5, col = "light gray", cex = 0.8)
  text(30, .1, paste("T50:", round(mean(na.omit(T50)), 1)), pos = 4, col = "gray")
  abline(v = round(mean(na.omit(T50)), 1), lty = 2, lwd = 1.5, col = "gray", cex = 0.8)
  text(30, .2, paste("T95:", round(mean(na.omit(T95)), 1)), pos = 4, col = "black")
  abline(v = round(mean(na.omit(T95)), 1), lty = 2, lwd = 1.5, col = "black", cex = 0.8)
}
