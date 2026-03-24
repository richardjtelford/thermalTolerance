#' Heat tolerances
#' @description 
#' Code for estimating heat tolerances using variable fluorescence!
#' @param temperature Temperature
#' @param fvfm fvfm data
#' @param id leaf ID
#' @param control_temp Control temperature
#' @param plot.est Should estimates be plotted
#' @param boots Number of bootstrap iterations
#' @examples
#' psiiht(temperature = htdata$temperature, fvfm = htdata$fvfm, 
#'   control_temp = 23, id = htdata$id, plot.est = TRUE, boots = 5)
#' @importFrom car logit
#' @importFrom stats nls coef lm na.omit quantile
#' @export


psiiht <- function(temperature, fvfm, control_temp, id, plot.est, boots) {
  HTdf <- data.frame(Temperature = temperature, fvfm = fvfm, id = id)

  return(do.call("rbind", by(HTdf, list(HTdf$id), function(df) {
    Temperature <- df[, which(colnames(df) == "Temperature")]
    fvfm <- df[, which(colnames(df) == "fvfm")]
    id <- df[, which(colnames(df) == "id")]
    # get parameter estimates for logistic decay model
    cof <- coef(lm(logit(fvfm) ~ Temperature))
    # Fit a non linear least squares model to the fvfm and Temperature data
    HT.model <- nls(fvfm ~ theta1 / (1 + exp(-(theta2 + theta3 * Temperature))),
      start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
      trace = FALSE, control = list(maxiter = 1000, tol = 1e-3)
    )

    # Use the parameter estimates (coef(HT.model)[#])from the HT.model to predict a new fit based on a heat treatments from 23-62 degrees celcius. Here, # = 1:3.
    y <- coef(HT.model)[1] / (1 + exp(-(coef(HT.model)[2] + coef(HT.model)[3] * seq(23, 62))))

    # Calculate half of the control Fv/Fm & a 95% reduction in fvfm with reference to control
    half <- mean(na.omit(fvfm[which(Temperature == control_temp)])) / 2
    nine5 <- mean(na.omit(fvfm[which(Temperature == control_temp)])) * 0.05
    # 95 Confidence Interval
    predict.boot <- matrix(NA, 40, boots)
    T95 <- T50 <- Tcrit <- c()
    for (k in 1:boots) {
      # print(k)
      srows <- sample(1:length(Temperature), length(Temperature), TRUE)

      if (inherits(try(nls(fvfm[srows] ~ theta1 / (1 + exp(-(theta2 + theta3 * Temperature[srows]))),
        start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
        trace = FALSE, control = list(maxiter = 1000, tol = 1e-3)
      ), silent = TRUE)[[1]], "nlsModel")) {
        HT.model2 <- nls(fvfm[srows] ~ theta1 / (1 + exp(-(theta2 + theta3 * Temperature[srows]))),
          start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
          trace = FALSE, control = list(maxiter = 1000, tol = 1e-6)
        )
        predict.boot[, k] <- coef(HT.model2)[1] / (1 + exp(-(coef(HT.model2)[2] + coef(HT.model2)[3] * seq(23, 62))))
        # Estimate T95
        T95[k] <- (-log((coef(HT.model2)[1] / nine5) - 1) - coef(HT.model2)[2]) / coef(HT.model2)[3]

        # Estimate T50
        T50[k] <- (-log((coef(HT.model2)[[1]] / half) - 1) - coef(HT.model2)[[2]]) / coef(HT.model2)[[3]]
        T50k <- (-log((coef(HT.model2)[[1]] / half) - 1) - coef(HT.model2)[[2]]) / coef(HT.model2)[[3]]
        # Use model to predict changes in fvfm & make new dataframe
        predict <- data.frame(x = seq(23, 62), y = coef(HT.model2)[1] / (1 + exp(-(coef(HT.model2)[2] + coef(HT.model2)[3] * seq(23, 62))))) # create a dataframe of predictions
        df1 <- cbind(predict[-1, ], predict[-nrow(predict), ])[, c(3, 1, 4, 2)]
        # Use new dataframe to estimate the slope at between each 1-degree interval
        df1$slp <- as.vector(apply(df1, 1, function(x) summary(lm((x[3:4]) ~ x[1:2]))[[4]][[2]]))
        slp.at.tcrit <- round(min(df1$slp), 3) * .15 # Determine where slope is 15% of max slope & round
        # Estimate the fvfm at which the slope is 15% of max slope & less than T50
        fvfv.at.tcrit <- df1[which(abs(df1[which(df1[, 1] < T50k), ]$slp - slp.at.tcrit) == min(abs(df1[which(df1[, 1] < T50k), ]$slp - slp.at.tcrit))), ][1, 3]
        Tcrit[k] <- (-log((coef(HT.model2)[[1]] / fvfv.at.tcrit) - 1) - coef(HT.model2)[[2]]) / coef(HT.model2)[[3]] # Estimate the temperatureat which the slope is 15% of max slope
      } else {
        (class(try(nls(fvfm ~ theta1 / (1 + exp(-(theta2 + theta3 * Temperature))),
          start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
          data = data2[srows, ], trace = FALSE, control = list(maxiter = 1000, tol = 1e-3)
        ), silent = TRUE)[[1]]) == "list")
        predict.boot[, k] <- NA
        T95[k] <- NA
        T50[k] <- NA
        Tcrit[k] <- NA
      }
    }

    fvfm.boot <- t(apply(predict.boot, 1, function(x) {
      quantile(x, c(0.025, 0.975), na.rm = TRUE)
    }))

    # Tcrit.ci=quantile(Tcrit,c(0.025,0.975),na.rm=TRUE)
    T50.ci <- quantile(T50, c(0.025, 0.975), na.rm = TRUE)
    # T95.ci=quantile(T95,c(0.025,0.975),na.rm=TRUE)

    if (plot.est) {
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
    return(data.frame(
      id = (unique(id)),
      # Tcrit.lci=round(Tcrit.ci[[1]],1),
      Tcrit.mn = round(mean(na.omit(Tcrit)), 1),
      # Tcrit.uci=round(Tcrit.ci[[2]],1),

      # T50.lci=round(T50.ci[[1]],1),
      T50.mn = round(mean(na.omit(T50)), 1),
      # T50.uci=round(T50.ci[[2]],1),

      # T95.lci=round(T95.ci[[1]],1),
      T95.mn = round(mean(na.omit(T95)), 1)
      # T95.uci=round(T95.ci[[2]],1)))
    ))
  })))
}
