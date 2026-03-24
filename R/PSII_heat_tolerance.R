#' Heat tolerances
#' @description
#' Code for estimating heat tolerances using variable fluorescence!
#' @param temperature Temperature
#' @param fvfm fvfm data
#' @param id leaf ID
#' @param control_temp Control temperature
#' @param boots Number of bootstrap iterations
#' @examples
#' psiiht(
#'   temperature = htdata$temperature, fvfm = htdata$fvfm,
#'   control_temp = 23, id = htdata$id, boots = 5
#' )
#' @importFrom car logit
#' @importFrom stats nls coef lm na.omit quantile
#' @importFrom purrr list_rbind 
#' @importFrom dplyr bind_cols
#' @export


psiiht <- function(temperature, fvfm, control_temp, id, boots) {
  HTdf <- data.frame(temperature = temperature, fvfm = fvfm, id = id) |>
    na.omit()
  
  # vector of temperatures to make predictions at
  t_vals <- seq(min(HTdf$temperature), max(HTdf$temperature), length = 100)

  split(HTdf, ~ id) |> 
    lapply(function(df) {
    # get parameter estimates for logistic decay model
    cof <- coef(lm(logit(fvfm) ~ temperature, data = df))
    # Fit a non linear least squares model to the fvfm and Temperature data
    HT.model <- nls(fvfm ~ theta1 / (1 + exp(-(theta2 + theta3 * temperature))),
      start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
      trace = FALSE, control = list(maxiter = 1000, tol = 1e-3),
      data = df
    )

    # Use the parameter estimates (coef(HT.model)[#])from the HT.model to predict a new fit based on a heat treatments from 23-62 degrees celcius. Here, # = 1:3.
    # y <- coef(HT.model)[1] / (1 + exp(-(coef(HT.model)[2] + coef(HT.model)[3] * seq(23, 62))))

    # Calculate half of the control Fv/Fm & a 95% reduction in fvfm with reference to control
    half <- mean(fvfm[which(df$temperature == control_temp)]) / 2
    nine5 <- mean(fvfm[which(df$temperature == control_temp)]) * 0.05
    # 95 Confidence Interval
    predict.boot <- matrix(NA, length(t_vals), boots)
    T95 <- T50 <- Tcrit <- numeric(boots)
    for (k in seq_len(boots)) {

      srows <- sample(x = 1:nrow(df), size = nrow(df), replace = TRUE)
      df_boot <- df[srows, ]
      HT.model2 <- try(
        nls(fvfm ~ theta1 / (1 + exp(-(theta2 + theta3 * temperature))),
          start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
          trace = FALSE, control = list(maxiter = 1000, tol = 1e-6),
          data = df_boot
        ),
        silent = TRUE
      )
      if (inherits(HT.model2[[1]], "nlsModel")) {
        predict.boot[, k] <- coef(HT.model2)[1] / (1 + exp(-(coef(HT.model2)[2] + coef(HT.model2)[3] * t_vals)))
        # Estimate T95
        T95[k] <- (-log((coef(HT.model2)[1] / nine5) - 1) - coef(HT.model2)[2]) / coef(HT.model2)[3]

        # Estimate T50
        T50[k] <- (-log((coef(HT.model2)[[1]] / half) - 1) - coef(HT.model2)[[2]]) / coef(HT.model2)[[3]]

        # Use model to predict changes in fvfm & make new dataframe
        predict <- data.frame(x = t_vals, y = coef(HT.model2)[1] / (1 + exp(-(coef(HT.model2)[2] + coef(HT.model2)[3] * t_vals)))) # create a dataframe of predictions
        df1 <- cbind(predict[-1, ], predict[-nrow(predict), ])[, c(3, 1, 4, 2)]
        # Use new dataframe to estimate the slope at between each 1-degree interval
        df1$slp <- as.vector(apply(df1, 1, function(x) summary(lm((x[3:4]) ~ x[1:2]))[[4]][[2]]))
        slp.at.tcrit <- round(min(df1$slp), 3) * .15 # Determine where slope is 15% of max slope & round
        # Estimate the fvfm at which the slope is 15% of max slope & less than T50
        fvfv.at.tcrit <- df1[which(abs(df1[which(df1[, 1] < T50[k]), ]$slp - slp.at.tcrit) == min(abs(df1[which(df1[, 1] < T50[k]), ]$slp - slp.at.tcrit))), ][1, 3]
        Tcrit[k] <- (-log((coef(HT.model2)[[1]] / fvfv.at.tcrit) - 1) - coef(HT.model2)[[2]]) / coef(HT.model2)[[3]] # Estimate the temperature at which the slope is 15% of max slope
      } else {
        predict.boot[, k] <- NA
        T95[k] <- NA
        T50[k] <- NA
        Tcrit[k] <- NA
      }
    }
    fvfm.boot <- apply(predict.boot, 1, function(x) {
      quantile(x, c(0.025, 0.975), na.rm = TRUE)
    }) |> 
      t() |> 
      as.data.frame() |> 
      set_names(c("low", "high")) |> 
      bind_cols(temperature = t_vals)

    list(
      data = df,
      fvfm.boot = fvfm.boot,
      T95 = T95,
      T50 = T50,
      Tcrit = Tcrit
    )
  })
}
