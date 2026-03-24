summary.htol <- function(htol) {

  
  # Tcrit.ci=quantile(Tcrit,c(0.025,0.975),na.rm=TRUE)
  T50.ci <- quantile(T50, c(0.025, 0.975), na.rm = TRUE)
  # T95.ci=quantile(T95,c(0.025,0.975),na.rm=TRUE)
  
    data.frame(
      id = unique(htol$df$id),
      # Tcrit.lci=round(Tcrit.ci[[1]],1),
      Tcrit.mn = round(mean(na.omit(htol$Tcrit)), 1),
      # Tcrit.uci=round(Tcrit.ci[[2]],1),

      # T50.lci=round(T50.ci[[1]],1),
      T50.mn = round(mean(na.omit(htol$T50)), 1),
      # T50.uci=round(T50.ci[[2]],1),

      # T95.lci=round(T95.ci[[1]],1),
      T95.mn = round(mean(na.omit(htol$T95)), 1)
      # T95.uci=round(T95.ci[[2]],1)))
    )
}