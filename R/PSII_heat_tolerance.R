#' Heat tolerances
#' @description
#' Code for estimating heat tolerances using variable fluorescence!
#' @param temperature Temperature
#' @param fvfm fvfm data
#' @param id leaf ID
#' @param control_temp Control temperature
#' @param boots Number of bootstrap iterations
#' @param warming Is the experiment analysing heat tolerance (TRUE)
#'  (fvfm declines with high values) or cold tolerance (FALSE)
#'  (fmfv declines with low values)? Important for calculating Tcrit.
#' @examples
#' library(future)
#' plan(multisession) # enable parallel processing
#' htol <- psiiht(
#'   temperature = htdata$temperature, fvfm = htdata$fvfm,
#'   control_temp = 23, id = htdata$id, boots = 5
#' )
#' summary(htol)
#' autoplot(htol)
#' @importFrom car logit
#' @importFrom stats nls coef lm na.omit quantile
#' @importFrom purrr list_rbind
#' @importFrom rlang set_names
#' @importFrom dplyr bind_cols filter slice_min pull
#' @importFrom futurize futurize
#' @export


psiiht <- function(temperature, fvfm, id, control_temp, warming = TRUE, boots) {
  HTdf <- data.frame(temperature = temperature, fvfm = fvfm, id = id) |>
    na.omit()

  # vector of temperatures to make predictions at
  t_vals <- seq(min(HTdf$temperature), max(HTdf$temperature), length = 100)

  results <- split(HTdf, ~id) |>
    lapply(function(df) {
      # get parameter estimates for logistic decay model
      cof <- coef(lm(logit(fvfm) ~ temperature, data = df))
      # Fit a non linear least squares model to the fvfm and Temperature data
      HT_model <- try(nls(fvfm ~ theta1 / (1 + exp(-(theta2 + theta3 * temperature))),
        start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
        trace = FALSE, control = list(maxiter = 1000, tol = 1e-3),
        data = df
      ))
      # return NULL if model does not fit
      if(inherits(HT_model, "try-error")) {
       return(NULL) 
      }

      # Use the parameter estimates (coef(HT_model)[#]) from the HT_model to
      # predict a new fit based on a heat treatments
      preds <- coef(HT_model)[1] /
        (1 + exp(-(coef(HT_model)[2] + coef(HT_model)[3] * t_vals)))

      # Calculate half of the control Fv/Fm & a 95% reduction in fvfm
      mean_control <- mean(df$fvfm[which(df$temperature == control_temp)])
      half <- mean_control / 2
      nine5 <- mean_control * 0.05
      # 95 Confidence Interval
      predict_boot <- matrix(NA, length(t_vals), boots)
      T95 <- T50 <- Tcrit <- numeric(boots)
      for (k in seq_len(boots)) {
        srows <- sample(x = nrow(df), size = nrow(df), replace = TRUE)
        df_boot <- df[srows, ]
        HT_model2 <- try(
          nls(fvfm ~ theta1 / (1 + exp(-(theta2 + theta3 * temperature))),
            start = list(theta1 = .8, theta2 = cof[1], theta3 = cof[2]),
            trace = FALSE, control = list(maxiter = 1000, tol = 1e-6),
            data = df_boot
          ),
          silent = TRUE
        )
        if (inherits(HT_model2[[1]], "nlsModel")) {
          predict_boot[, k] <- coef(HT_model2)[1] /
            (1 + exp(-(coef(HT_model2)[2] + coef(HT_model2)[3] * t_vals)))
          # Estimate T95
          T95[k] <- (-log((coef(HT_model2)[1] / nine5) - 1) -
            coef(HT_model2)[2]) / coef(HT_model2)[3]

          # Estimate T50
          T50[k] <- (-log((coef(HT_model2)[[1]] / half) - 1) -
            coef(HT_model2)[[2]]) / coef(HT_model2)[[3]]

          # Use model to predict changes in fvfm & make new dataframe
          # create a dataframe of predictions

          predict <- data.frame(
            x = t_vals,
            y = coef(HT_model2)[1] /
              (1 + exp(-(coef(HT_model2)[2] + coef(HT_model2)[3] * t_vals)))
          )
          df1 <- cbind(predict[-1, ], predict[-nrow(predict), ])[, c(3, 1, 4, 2)]
          # Use new dataframe to estimate the slope at between each interval
          df1$slp <- as.vector(apply(df1, 1, function(x) {
            coef(lm((x[3:4]) ~ x[1:2]))[2]
          }))
          # Determine where slope is 15% of max slope & round
          if(warming) {
            slp.at.tcrit <- min(df1$slp) * .15
            # Find the temperature at which the slope is 15% of max slope and < T50
            Tcrit[k] <- df1 |> 
              filter(.data$x < T50[k]) |> 
              slice_min(abs(.data$slp - slp.at.tcrit)) |> 
              pull(.data$x)
          } else { # cooling
            slp.at.tcrit <- max(df1$slp) * .15
            # Find the temperature at which the slope is 15% of max slope and > T50
            
            Tcrit[k] <- df1 |> filter(.data$x > T50[k]) |> 
              slice_min(abs(.data$slp - slp.at.tcrit)) |> 
              pull(.data$x)
            
          }
        } else {
          predict_boot[, k] <- NA
          T95[k] <- NA
          T50[k] <- NA
          Tcrit[k] <- NA
        }
      }
      preds <- apply(predict_boot, 1, function(x) {
        quantile(x, c(0.025, 0.975), na.rm = TRUE)
      }) |>
        t() |>
        as.data.frame() |>
        set_names(c("low", "high")) |>
        bind_cols(temperature = t_vals, preds = preds)

      list(
        preds = preds,
        T95 = T95,
        T50 = T50,
        Tcrit = Tcrit
      )
    }) |> futurize(seed = TRUE)
  results <- list(data = HTdf, results = results)
  class(results) <- c("htol", class(results))
  results
}
