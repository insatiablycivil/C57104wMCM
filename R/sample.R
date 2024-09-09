#' Sample: Randomly generate samples to later apply Screening to for Monte Carlo estimation
#'
#' This function generates samples for the Monte Carlo based on the statistical description provided.
#'
#' @param gas_data Data.frame: First column containing days since oldest sample. Rest of columns containing sample value in PPM.
#' @param limit_data Data.frame: Columns containing the IEEE DGA Table limits for each gas. Each Table is on its own row; four rows expected. Caution: Table 4 is expressed as PPM/Day!
#' @param accuracy Named Numeric: Values to be multiplied by gas value to give uncertainty. Can either be same value for all gases, or individual values for each gas. Expected range: 0-1.
#' @param min_accuracy Named Numeric: Absolute Values (PPM) used as minimum uncertainty. Can either be same value for all gases, or individual values for each gas.
#' @param distribution_shape Numeric: Distribution shape representing uncertainty. 1 = Normal, 2 = Triangular, 3 = Uniform. If Normal, Accuracy is interpreted as 95% confidence interval, else as 100%.
#' @param rho Numeric: Value representing rho for correlations. Can either be same value for all gases, or individual values for each gas. NA if using \code{cor_mat} instead.
#' @param cor_mat Numeric Matrix: Values representing correlation matrix. NA if using \code{rho} instead.
#' @param n Numeric: Number of samples to generate for the Monte Carlo for each gas.
#' @param col_names Strings: Gas Names. Order must match other named variables throughout. Expected a subset of <H2|CH4|C2H6|C2H4|C2H2|CO|CO2>.
#' @param sd_data Data.frame: First column containing days since oldest sample. Rest of columns containing standard deviation for sample value in PPM.
#' @return Numeric Array: Values representing generated samples for the Monte Carlo based on the statistical description provided. Shaped: Days x Gases x Samples.
#' @examples
#' Sample(gas_data = data.frame(date = c(4,3,2,1,0), H2 = c(8,4,2,1,0.5), CH4 = c(10,20,40,80,160)),
#'        limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
#'        accuracy = setNames(c(0.15, 0.50), c("H2", "CH4")),
#'        min_accuracy = setNames(c(1, 5), c("H2", "CH4")),
#'        distribution_shape = 1,
#'        rho = 0.5,
#'        cor_mat = NA,
#'        n = 1000,
#'        col_names = c("H2", "CH4"),
#'        sd_data = data.frame(date = c(4,3,2,1,0), H2 = c(1.2,0.6,0.3,0.15,0.075), CH4 = c(5,10,20,40,80)))
#' Sample(gas_data = data.frame(date = c(4,3,2,1,0), H2 = c(8,4,2,1,0.5), CH4 = c(10,20,40,80,160)),
#'        limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
#'        accuracy = 0.15,
#'        min_accuracy = 1,
#'        distribution_shape = 1,
#'        rho = NA,
#'        cor_mat = matrix(c(1,0,0,1), ncol=2),
#'        n = 1000,
#'        col_names = c("H2", "CH4"),
#'        sd_data = data.frame(date = c(4,3,2,1,0), H2 = c(1.2,0.6,0.3,0.15,0.075), CH4 = c(0.75,1.50,3,6,12)))
#' @export
Sample <- function(gas_data, limit_data, accuracy, min_accuracy,
                   distribution_shape, rho, cor_mat, n, col_names, sd_data) {
  # Generate samples
  d = nrow(gas_data) # number of samples
  g = length(col_names) # number of gases

  # Need permutations for different stuff
  mus <- unlist(gas_data[,col_names])
  mus_g_d <- as.numeric(t(array(mus, dim = c(d, g)))) # val for each gas for each day
  sds <- unlist(sd_data[,col_names])
  sd_g_d  <- as.numeric(t(array(sds, dim = c(d, g)))) # SD for each gas for each day

  if((all(!is.na(rho)) && (sum(rho)==0)) |
     (all(!is.na(cor_mat)) && (sum(cor_mat)==sum(diag(cor_mat))))) { # Independent
    if(distribution_shape == 1) { # Normal Dist, assumed Accuracy given to 95%
      y_N <- rnorm(n*d*g, mus, sds)
    } else {
      # Assumes Accuracy given as limits (100%)
      # Fiddly reordering stuff from grouped by gas to grouped by sample (then back)
      lims <- pmax(mus_g_d*accuracy, min_accuracy)
      lims <- as.numeric(t(array(lims, dim= c(g,d))))
      y_N <- runif(n*d*g, mus-lims, mus+lims) # Uniform
      if(distribution_shape == 2) { # Triangle
        mask <- (as.numeric(y_N < mus)*2) - 1 # LHS == 1, RHS == -1
        y_N <- mus + mask*(sqrt(abs(y_N-mus)*lims)-lims)
      }
    }
    y_N <- array(y_N, dim = c(d, g, n)) # days x gases x samples
  } else { # Dependence
    library(MASS) # needed for mvrnorm
    if(all(!is.na(rho))) {
      rho_mat <- (matrix((rho), nrow = g, ncol = g) + diag(rep((1-rho),g)))
      rho_mat <- diag(1, d) %x% rho_mat
    } else {
      rho_mat <- diag(1,d) %x% cor_mat
    }
    cor_mat <- diag(sd_g_d)
    sigma <- cor_mat %*% rho_mat %*% cor_mat
    # Idea from (stats.stackexchange.com/q/143280): Start with Norm, then convert
    y_N <- array(mvrnorm(n, mus_g_d, sigma), dim = c(n, g, d))
    y_N <- aperm(y_N, c(2, 3, 1)) # gases x days x samples
    if(distribution_shape != 1) { # Need to convert distribution shape
      y_N <- pnorm(as.numeric(y_N), mus_g_d, sd_g_d) # Uniform
      if(distribution_shape == 2) { # Triangle
        mask <- (as.numeric(y_N < 0.5)*2) - 1 # LHS == 1, RHS == -1
        y_N <- 0.5 + mask*(sqrt(abs(y_N-0.5)*0.5)-0.5)
      }
      # Assumes Accuracy given as limits (100%)
      lims <- pmax(mus_g_d*accuracy, min_accuracy)
      y_N <- (((y_N-0.5)*2)*lims) + mus_g_d # Rescaling after shape change
      y_N <- array(y_N, dim = c(g,d,n)) # gases x days x samples
    }
    y_N <- aperm(y_N, c(2, 1, 3)) # days x gases x samples
  }
  return(y_N)
}
