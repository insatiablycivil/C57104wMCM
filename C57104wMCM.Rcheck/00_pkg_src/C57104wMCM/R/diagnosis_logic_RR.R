#' Interpret_RR: Applies Rogers Ratio Diagnostic Logic as per IEEE Table 5
#'
#' This functions interprets the inputted gas sample ratios to determine the appropriate diagnosis for each in accordance to IEEE C57.104-2019 Table 5.
#'
#' @param y_in Numeric Array: Values representing generated gas ratios, assumed order: C2H2/C2H4, CH4/H2, C2H4/C2H6, assumed shape: 3 x Samples.
#' @param n_in Numeric: Number of samples to be diagnosed.
#' @return Strings: List of Diagnoses, one for each inputted sample.
#' @examples
#' Interpret_RR(y_in = array(c(0.090, 0.002, 0.401, 0.120, 0.001, 0.310), dim = c(3,2)),
#'              n_in = 2)
#' @export
Interpret_RR <- function(y_in, n_in) {
  # Assumed y_n order: C2H2/C2H4, CH4/H2, C2H4/C2H6
  # As in C57.104-2019 Table 5.
  if(n_in == 1) y_in <- array(rep(unlist(y_in), 2), dim=c(3,2)) # Add dummy column for shape
  D = rep("ND", n_in)
  #         Lower Limit          Upper Limit
  i_D =                         (y_in[1,] <  0.1) &
           (y_in[2,] >= 0.1) &  (y_in[2,] <= 1.0) &
                                (y_in[3,] <  1.0)
  D[i_D] = "C0"
  i_D =                         (y_in[1,] <  0.1) &
                                (y_in[2,] <  0.1) &
                                (y_in[3,] <  1.0)
  D[i_D] = "C1"
  i_D =    (y_in[1,] >= 0.1) &  (y_in[1,] <= 3.0) &
           (y_in[2,] >= 0.1) &  (y_in[2,] <= 1.0) &
           (y_in[3,] >  3.0)
  D[i_D] = "C2"
  i_D =                         (y_in[1,] <  0.1) &
           (y_in[2,] >= 0.1) &  (y_in[2,] <= 1.0) &
           (y_in[3,] >= 1.0) &  (y_in[3,] <= 3.0)
  D[i_D] = "C3"
  i_D =                         (y_in[1,] <  0.1) &
           (y_in[2,] >  1.0) &
           (y_in[3,] >= 1.0) &  (y_in[3,] <= 3.0)
  D[i_D] = "C4"
  i_D =                         (y_in[1,] <  0.1) &
           (y_in[2,] >  1.0) &
           (y_in[3,] >  3.0)
  D[i_D] = "C5"

  if(n_in == 1) return(D[1])
  return(D)
}
