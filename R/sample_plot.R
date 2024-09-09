#' Sample_Plot: Plotting outputs from Sample function for validation.
#'
#' This is internal function to plot generated samples for validation.
#'
#' @param y_N Numeric Array: Values representing generated samples for the Monte Carlo based on the statistical description provided. Shaped: Days x Gases x Samples.
#' @param PP List: Outputs from Preprocess function.
#' @examples
#' # Not recreating required inputs to run this example
#' # Sample_Plot(y_1 = y_N[1,,],
#' #             PP = preprocessed)
#' @export
Sample_Plot <- function(y_N, PP) {
  # For viewing the data for validation purposes.
  # y_N has samples [days,gases,samples]
  # PP (preprocessed) has various parameters
  # Looking at:
  #   1: gas values over time to check input is as expected
  #   2: gas values over time with uncertainty to show how that was interpreted
  #   3: correlation plot of gas samples for each day to show:
  #       a: marginal distribution shape is as expected
  #       b: correlations across gases are as expected

  library(tidyverse)
  library(GGally)

  Plot_Uni_Dist <- function(x, title, x_label, y) {
    d_x <- seq(min(x), max(x), length.out = bin_count + 1)
    mu_in <- df$gas_value[y]
    sd_in <- df$sd_value[y]
    lim_in <- df$lim[y]
    w <- ((max(x)-min(x)) / bin_count)*n
    hist(x, main = title, xlab = x_label, breaks = d_x)
    if(PP$distribution_shape == 1) d_density <- dnorm(d_x, mean=mu_in, sd=sd_in)*w
    if(PP$distribution_shape == 2) d_density <- dtriangle(d_x, a=mu_in-lim_in, b=mu_in+lim_in)*w
    if(PP$distribution_shape == 3) d_density <- dunif(d_x, min=mu_in-lim_in, max=mu_in+lim_in)*w
    lines(d_x, d_density, col = "red")
  }

  Diag_Plot_Helper <- function(data, mapping, ...) {
    x <- unlist(lapply(mapping, all.vars))
    i <- match(x, PP$col_names)
    w <- ((max(data[,i])-min(data[,i])) / bin_count)*n
    mask <- (df$date == PP$gas_data$date[idx]) & (df$Gas == x)
    mu_in <- df$gas_value[mask]
    sd_in <- df$sd_value[mask]
    lim_in <- df$lim[mask]
    p <- ggplot(data = data, mapping = mapping) +
         geom_histogram(bins = bin_count)
    if(PP$distribution_shape == 1) p <- p+stat_function(fun=function(x){dnorm(x, mean=mu_in, sd=sd_in)*w}, colour="red")
    if(PP$distribution_shape == 2) p <- p+stat_function(fun=function(x){dtriangle(x, a=mu_in-lim_in, b=mu_in+lim_in)*w}, colour="red")
    if(PP$distribution_shape == 3) p <- p+stat_function(fun=function(x){dunif(x, min=mu_in-lim_in, max=mu_in+lim_in)*w}, colour="red")
    p
  }

  Lower_Plot_Helper <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(colour = "darkgrey", alpha = 0.5) +
      geom_smooth(method = "lm", color = "red", ...)
    p
  }

  d <- dim(y_N)[1]
  g <- dim(y_N)[2]
  n <- dim(y_N)[3]

  # Get Gas Values and Time
  df <- PP$gas_data %>% pivot_longer(cols = -date, names_to = "Gas", values_to = "gas_value")
  df$COn <- c("CnHn","COn")[as.numeric(df$Gas %in% c("CO","CO2"))+1] # Seperating CnHn and COn
  # Plotting 1: Gas Samples Over Time
  suppressWarnings(suppressMessages(print(ggplot(data = df, aes(x = date, y = gas_value, color = Gas)) +
        geom_line() +
        geom_point(aes(shape = Gas)) +
        labs(x = "Days Since Oldest Sample", y = "Gas Value [PPM]") +
        facet_wrap(~ COn, ncol = 1, scales = "free_y") +
        theme_minimal() + theme(legend.position="bottom") +
        ggtitle("Input Gas Sample Values"))))

  # Get SD, calculate limit to see if min_acc was used.
  df <- PP$sd_data %>% pivot_longer(cols = -date, names_to = "Gas", values_to = "sd_value") %>%
    left_join(df, by = c("date", "Gas"))
  sub_title <- ifelse(PP$distribution_shape == 1,
                     "Normal Distribution: Bounds showing 95% confidence",
                     paste0(c("Triangle","Uniform")[PP$distribution_shape], " Distribution: Bounds showing absolute limits"))
  df$lim_u <- df$gas_value*PP$accuracy
  if(length(PP$min_accuracy) == 1) {
    df$lim_u_min <- PP$min_accuracy
  } else {
    df$lim_u_min <- rep(PP$min_accuracy, d)
  }
  df$lim <- pmax(df$gas_value*PP$accuracy, PP$min_accuracy)

  # Plotting 2: Gas Samples with Uncertainty Over Time (If Normal, use 95% else treat as limit)
  suppressWarnings(suppressMessages(print(ggplot(data = df, aes(x = date)) +
      geom_line(aes(y = gas_value, color = Gas)) +
      geom_errorbar(data = df %>% filter(lim==lim_u),
                    aes(ymin = gas_value-lim, ymax = gas_value+lim), color = "black", width=.1) +
      geom_errorbar(data = df %>% filter(lim==lim_u_min),
                    aes(ymin = gas_value-lim, ymax = gas_value+lim), color = "red", width=.1) +
      labs(x = "Days Since Oldest Sample", y = "Gas Value [PPM]") +
      facet_wrap(~ Gas, ncol = 1, scales = "free_y") +
      theme_minimal() + theme(legend.position="none") +
      ggtitle("Input Gas Sample Values with Uncertainty",
              subtitle = paste0(sub_title, "\nRed Error Bars Indicate Minimum Accuracy Overrode")))))

  # Plotting 3: Correlations across Gas Samples per Day
  bin_count <- ceiling(log(n)*2)
  for(idx in c(1:d)) {
    if(g == 1) {
      Plot_Uni_Dist(t(y_N[idx,,]),
            paste0("Distribution of Generated Gas Samples: Day ", PP$gas_data$date[idx]),
            paste0(PP$col_names, " [ppm]"),
            idx)
    } else {
      df2 <- as.data.frame(t(y_N[idx,,]))
      colnames(df2) <- PP$col_names
      suppressWarnings(suppressMessages(ggpairs(df2,
                    diag = list(continuous = Diag_Plot_Helper),
                    lower = list(continuous = Lower_Plot_Helper)) +
            ggtitle(paste0("Correlations Across Generated Gas Samples: Day ", PP$gas_data$date[idx]))))
    }
  }
}
