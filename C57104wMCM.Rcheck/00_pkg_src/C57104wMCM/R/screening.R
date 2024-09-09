#' Screening: Apply Simplified IEEE C57.104-2019 Screening - Use all samples for Tables 1-4.
#'
#' This function applies simplified IEEE C57.104-2019 Screening. It uses all samples for Tables 1-4 to give a per-gas and combined DGA Status. The frequency of given outcomes is interpreted as probabilities of them occurring.
#'
#' @param y_N Numeric Array: Values representing generated samples for the Monte Carlo based on the statistical description provided. Shaped: Days x Gases x Samples.
#' @param gas_data Data.frame: First column containing days since oldest sample. Rest of columns containing sample value in PPM.
#' @param limit_data Data.frame: Columns containing the IEEE DGA Table limits for each gas. Each Table is on its own row; four rows expected. Caution: Table 4 is expressed as PPM/Day!
#' @param silent Boolean: Defaulted to FALSE.
#' @param verbose Boolean: Defaulted to FALSE.
#' @return List containing:
#' \itemize{
#'  \item $\code{protocol}: Numeric. 1|2|3 depending on if 1, 2, or more samples.
#'  \item $\code{T<n>_default}: Output for Table <n> using default gas values.
#'  \item $\code{T<n>}: Probability of output for Table <n> using generated samples.
#'  \item $\code{T<n>_bounds}: Numeric. Associated estimate of confidence interval.
#'  \item $\code{gas_L<n>_default}: Numeric. Output for per-gas DGA Status Level <n>.
#'  \item $\code{gas_L<n>}: Numeric. Probability of output for per-gas DGA Status Level <n>.
#'  \item $\code{gas_L<n>_bounds}: Numeric. Associated confidence interval.
#'  \item $\code{L<n>_default}: Numeric. Output for combined DGA Status Level <n>.
#'  \item $\code{L<n>}: Numeric. Probability of output for combined DGA Status Level <n>.
#'  \item $\code{L<n>_bounds}: Numeric. Mask indicated which generated samples were combined DGA Status Level <n>.
#'  \item $\code{L<n>_mask}: Logical Array. Associated estimate of confidence interval.
#'  \item $\code{T12_vals_default}: Numeric. Metric values using default gas values used for Tables 1 and 2.
#'  \item $\code{T3_vals_default}: Numeric. Metric values using default gas values used for Table 3.
#'  \item $\code{T4_vals_default}: Numeric. Metric values using default gas values used for Tables 4.
#'  \item $\code{T12_vals}: Numeric. Metric values using generated samples used for Tables 1 and 2.
#'  \item $\code{T3_vals}: Numeric. Metric values using generated samples used for Table 3.
#'  \item $\code{T4_vals}: Numeric. Metric values using generated samples used for Tables 4.
#'  \item $\code{GL<n>_mask}: Logical Array. Mask indicated which generated samples were per-gas DGA Status Level <n>.
#' }
#' @examples
#' Screening(y_N = array(c(2,1,0.5,40,80,160,2.1,1.1,0.6,40.1,80.1,160.1), dim = c(2,2,2)),
#'           gas_data = data.frame(date = as.Date(c(2,1,0)), H2 = c(8,4,1), CH4 = c(40,80,160)),
#'           limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
#'           silent = FALSE, verbose = FALSE)
#' @export
Screening <- function(y_N, limit_data, gas_data, silent = FALSE, verbose = FALSE) {
  # Applies simplified IEEE Screening
  # Three protocols: 1 -> 1 sample, 2 -> 2 samples, 3 -> >3 samples
  library(abind)

  bounds <- function(p) {
    return(((1.96 * sqrt((p*(1-p))/n)))*100)
  }

  if(silent == FALSE) {
    print_message <- function(label, value, opt_val = NULL) {
      head <- sprintf("%8s", label)
      body <- sprintf("%07.3f", value)
      cat(head, ":<", body, ">\n")
      if(!is.null(opt_val)) {
        head <- sprintf("%8s", "+/-")
        body <- sprintf("%07.3f", opt_val)
        cat(head, ":<", body, ">\n")
      }
    }
  }

  X_list <- cbind(gas_data$date) # X axis (dates) won't be changing
  X_mat  <- cbind(1, X_list)
  d <- dim(y_N)[1] # y_N is [days, gas, samples]
  g <- dim(y_N)[2]
  n <- dim(y_N)[3]
  T1 <- as.numeric(limit_data[1,, drop = FALSE])
  T2 <- as.numeric(limit_data[2,, drop = FALSE])
  if(d > 1) T3 <- as.numeric(limit_data[3,, drop = FALSE])
  if(d > 2) T4 <- as.numeric(limit_data[4,, drop = FALSE])

  # Tables 1|2:
  #   Getting Defaults
  D_V_T12 = gas_data[1,-1] # This is a bit redundant
  D_T1 = D_V_T12 < T1
  D_T2 = D_V_T12 < T2
  #   Getting Pass Indices
  V_T12 = y_N[1,,, drop = FALSE] # This is a bit redundant
  i_T1 = (V_T12 < T1)
  dim(i_T1) <- c(g, n)
  i_T2 = (V_T12 < T2)
  dim(i_T2) <- c(g, n)
  #   Getting Counts
  n_T1 = rowSums(i_T1)
  n_T2 = rowSums(i_T2)
  n_T2n1 = rowSums((!i_T1) & i_T2)
  #   Getting Probabilities
  P_T1 = n_T1/n
  P_T2 = n_T2/n
  P_T2n1 = n_T2n1/n
  #   Getting Estimated Uncertainty Bounds
  U_T1 = bounds(P_T1)
  U_T2 = bounds(P_T2)
  U_T2n1 = bounds(P_T2n1)
  # Tables 1|2|3:
  if(d > 1) {
    #   Getting Default
    D_V_T3 = (gas_data[1,-1] - gas_data[2,-1])
    D_T3 = D_V_T3 < T3
    #   Getting Metric
    V_T3 = y_N[1,,, drop = FALSE] - y_N[2,,, drop = FALSE] # [gases, samples]
    dim(V_T3) <- c(g, n)
    #   Getting Pass Indices
    i_T3 = (V_T3 < T3)
    dim(i_T3) <- c(g, n)
    #   Getting Counts
    n_T3 = rowSums(i_T3)
    n_T1n3 = ifelse((n_T1==0)|(n_T3==n), rep(0,g), rowSums(i_T1 & (!i_T3)))
    #   Getting Probabilities
    P_T3 = n_T3/n
    P_T1n3 = n_T1n3/n
    #   Getting Estimated Uncertainty Bounds
    U_T3 = bounds(P_T3)
    U_T1n3 = bounds(P_T1n3)
  } else {
    D_V_T3 <- NA
    D_T3 <- NA
    V_T3 <- NA
    P_T3 <- NA
    U_T3 <- NA
  }
  # Tables 1|2|3|4:
  if(d > 2) {
    #   Getting Default
    D_V_T4 <- lm.fit(X_mat, array(unlist(gas_data[,-1]), dim = c(d,g)))$coefficients
    if(g == 1) { # Output shape depends on whether 1 gas or not
      D_V_T4 <- D_V_T4[2]
    } else {
      D_V_T4 <- D_V_T4[2,]
    }
    D_T4 = D_V_T4 < T4
    #   Getting Metric
    # lm.fit is meant for 2D. Going to stack gasses
    y_N = aperm(y_N, c(1, 3, 2))
    y_N = array(y_N, dim = c(d, g*n))
    V_T4 <- lm.fit(X_mat, y_N)$coefficients[2,]
    V_T4 <- array(V_T4, dim = c(n, g))
    V_T4 = aperm(V_T4, c(2,1))
    y_N <- array(y_N, dim=c(d,n,g))
    y_N = aperm(y_N, c(1, 3, 2))
    #   Getting Pass Indices
    i_T4 = (V_T4 < T4)
    #   Getting Counts (Some are redundant here)
    n_T4 = rowSums(i_T4)
    n_Tn4 = rowSums(!i_T4)
    n_T24n1 = ifelse((n_T1==n)|(n_T2==0)|(n_T4==0), 0, rowSums((!i_T1) & i_T2 & i_T4))
    n_T4n2 = ifelse((n_T2==n)|(n_T4==0), 0, rowSums((!i_T2) & i_T4))
    n_T234n1 = ifelse((n_T1==n)|(n_T2==0)|(n_T3==0)|(n_T4==0), 0, rowSums((!i_T1) & i_T2 & i_T3 & i_T4))
    n_T24n3 = ifelse((n_T2==0)|(n_T3==n)|(n_T4==0), 0, rowSums(i_T2 & (!i_T3) & i_T4))
    #   Getting Probabilities
    P_T4 = n_T4/n
    P_Tn4 = n_Tn4/n
    P_T24n1 = n_T24n1/n
    P_T4n2 = n_T4n2/n
    P_T234n1 = n_T234n1/n
    P_T24n3 = n_T24n3/n
    #   Getting Estimated Uncertainty Bounds
    U_T4 = bounds(P_T4)
    U_Tn4 = bounds(P_Tn4)
    U_T24n1 = bounds(P_T24n1)
    U_T4n2 = bounds(P_T4n2)
    U_T234n1 = bounds(P_T234n1)
    U_T24n3 = bounds(P_T24n3)
  } else {
    D_V_T4 <- NA
    D_T4 <- NA
    V_T4 <- NA
    P_T4 <- NA
    U_T4 <- NA
  }
  # Per Gas Status Level:
  #   Getting Default and Indices
  if(d == 1) { # Protocol 1: Tables 1|2:
    i_GL3 = (!i_T2)
    i_GL1 = (i_T1)
    i_GL2 = (!i_GL3) & (!i_GL1)
    D_GL3 = (!D_T2)
    D_GL1 = (D_T1)
    D_GL2 = (!D_GL3) & (!D_GL1)
  }
  if(d == 2) { # Protocol 2: Tables 1|2|3:
    i_GL3 = (!i_T2)
    i_GL1 = (i_T1 & i_T3)
    i_GL2 = (!i_GL3) & (!i_GL1)
    D_GL3 = (!D_T2)
    D_GL1 = (D_T1 & D_T3)
    D_GL2 = (!D_GL3) & (!D_GL1)
  }
  if(d > 2) { # Protocol 3: Tables 1|2|3|4:
    i_GL3 = (i_T4 & (!i_T2)) | (!i_T4)
    i_GL1 = (i_T1 & i_T2 & i_T3 & i_T4)
    i_GL2 = (!i_GL3) & (!i_GL1)
    D_GL3 = (D_T4 & (!D_T2)) | (!D_T4)
    D_GL1 = (D_T1 & D_T2 & D_T3 & D_T4)
    D_GL2 = (!D_GL3) & (!D_GL1)
  }
  #   Getting Probabilities
  P_GL1 = rowSums(i_GL1)/n
  P_GL2 = rowSums(i_GL2)/n
  P_GL3 = rowSums(i_GL3)/n
  #   Getting Estimated Uncertainty Bounds
  U_GL1 = bounds(P_GL1)
  U_GL2 = bounds(P_GL2)
  U_GL3 = bounds(P_GL3)
  # Combined Status Level:
  #   Getting Default and Indices
  if(g == 1) {
    i_SL1 = i_GL1
    i_SL2 = i_GL2
    i_SL3 = i_GL3
    D_SL1 = D_GL1
    D_SL2 = D_GL2
    D_SL3 = D_GL3
  } else { # Pick Worst Case
    i_SL3 = as.logical(ceiling(colMeans(i_GL3)))
    i_SL1 = as.logical(floor(colMeans(i_GL1)))
    i_SL2 = (!i_SL3) & (!i_SL1)
    D_SL3 = as.logical(max(D_GL3))
    D_SL1 = as.logical(min(D_GL1))
    D_SL2 = (!D_SL3) & (!D_SL1)
  }
  #   Getting Probabilities
  P_SL1 = sum(i_SL1)/n
  P_SL2 = sum(i_SL2)/n
  P_SL3 = sum(i_SL3)/n
  #   Getting Estimated Uncertainty Bounds
  U_SL1 = bounds(P_SL1)
  U_SL2 = bounds(P_SL2)
  U_SL3 = bounds(P_SL3)
  #   Printing Outputs
  if(silent == FALSE) {
    cat("\nDefault Output: Combined Status Level: ",
        c("SL1","SL2","SL3")[D_SL1+(2*D_SL2)+(3*D_SL3)], "\n")
    cat("\nAll Samples: Combined Status Levels:\n")
    print_message("P SL1", P_SL1*100, U_SL1)
    print_message("P SL2", P_SL2*100, U_SL2)
    print_message("P SL3", P_SL3*100, U_SL3)
    if(verbose == TRUE) {
      cat("\nDefault Output: Per Gas Tables:\n")
      cat("Table 1: <", c("PASS","FAIL")[as.numeric(D_T1==0)+1], ">\n")
      cat("Table 2: <", c("PASS","FAIL")[as.numeric(D_T2==0)+1], ">\n")
      if(d > 1) cat("Table 3: <", c("PASS","FAIL")[as.numeric(D_T3==0)+1], ">\n")
      if(d > 2) cat("Table 4: <", c("PASS","FAIL")[as.numeric(D_T4==0)+1], ">\n")
      cat("\nAll Samples: Per Gas Tables:\n")
      print_message("P T1", P_T1*100, U_T1)
      print_message("P T2", P_T2*100, U_T2)
      if(d > 1) print_message("P T3", P_T3*100, U_T3)
      if(d > 2) print_message("P T4", P_T4*100, U_T4)
      if(g == 1) {
        print_message("P T2n1", P_T2n1*100, U_T2n1)
      }
      if(g == 2) {
        print_message("P T1n3", P_T1n3*100, U_T1n3)
      }
      if(g > 2) {
        print_message("P T234n1", P_T234n1*100, U_T234n1)
        print_message("P T24n3", P_T24n3*100, U_T24n3)
        print_message("P T4n2", P_T4n2*100, U_T4n2)
        print_message("P Tn4", P_Tn4*100, U_Tn4)
      }
    }
    if(g > 1) {
      cat("\nDefault Output: Per Gas Status Level: <",
          c("GL1","GL2","GL3")[D_SL1+(2*D_SL2)+(3*D_SL3)], ">\n")
      cat("\nper gas Status Levels:\n")
      print_message("P GL1", P_GL1*100, U_GL1)
      print_message("P GL2", P_GL2*100, U_GL2)
      print_message("P GL3", P_GL3*100, U_GL3)
    }
  }

  output <- list(protocol = which(c(g==1,g==2,g>2)),
                 T1_default = D_T1,
                 T2_default = D_T2,
                 T3_default = D_T3,
                 T4_default = D_T4,
                 T1 = P_T1,
                 T2 = P_T2,
                 T3 = P_T3,
                 T4 = P_T4,
                 T1_bounds = U_T1,
                 T2_bounds = U_T2,
                 T3_bounds = U_T3,
                 T4_bounds = U_T4,
                 gas_L1_default = D_GL1,
                 gas_L2_default = D_GL2,
                 gas_L3_default = D_GL3,
                 gas_L1 = P_GL1,
                 gas_L2 = P_GL2,
                 gas_L3 = P_GL3,
                 gas_L1_bounds = U_GL1,
                 gas_L2_bounds = U_GL2,
                 gas_L3_bounds = U_GL3,
                 L1_default = D_SL1,
                 L2_default = D_SL2,
                 L3_default = D_SL3,
                 L1 = P_SL1,
                 L2 = P_SL2,
                 L3 = P_SL3,
                 L1_bounds = U_SL1,
                 L2_bounds = U_SL2,
                 L3_bounds = U_SL3,
                 L1_mask = i_SL1,
                 L2_mask = i_SL2,
                 L3_mask = i_SL3,
                 T12_vals_default = D_V_T12,
                 T3_vals_default = D_V_T3,
                 T4_vals_default = D_V_T4,
                 T12_vals = V_T12,
                 T3_vals = V_T3,
                 T4_vals = V_T4,
                 GL1_mask = i_GL1,
                 GL2_mask = i_GL2,
                 GL3_mask = i_GL3)
  return(output)
}
