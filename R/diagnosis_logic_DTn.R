#' Interpret_DTn: Applies Duval Triangles 1, 4, or 5 Diagnostic Logic as per IEEE Table 6, D.3, and D.4, respectively.
#'
#' This functions interprets the inputted gas sample ratios to determine the appropriate diagnosis for each in accordance to either IEEE C57.104-2019 Table 6 or Table D.3, or Table D.4. Table D.3 is modified here such that original definition for C: (CH4 >= 36 & C2H6 >= 24) is instead: (CH4 >= 36 & C2H6 < 24).
#'
#' @param DT_n String: Value corresponding to either Triangle 1, 4, or 5. Expected value: 1, 4, or 5.
#' @param y_in Numeric Array: Values representing generated gas ratios, assumed order: DT1: C2H4, CH4, C2H2, DT4: CH4, H2, C2H6, DT5: C2H4, CH4, C2H6. Assumed shape: 3 x Samples.
#' @param n_in Numeric: Number of samples to be diagnosed.
#' @return Strings: List of Diagnoses, one for each inputted sample.
#' @examples
#' Interpret_DTn(DT_n = "4",
#'                y_in = array(c(0.090, 0.002, 0.401, 0.120, 0.001, 0.310), dim = c(3,2)),
#'                n_in = 2)
Interpret_DTn <- function(DT_n, y_in, n_in) {
  # Assumed y_n order (DIFFERS TO TABLE ORDER):
    # DT1: C2H4, CH4, C2H2
    # DT4: CH4, H2, C2H6
    # DT5: C2H4, CH4, C2H6
  Interpret_DT1 <- function(y_in, n_in) {
    # As in C57.104-2019 Table 6
    if(n_in == 1) y_in <- array(rep(unlist(y_in), 2), dim=c(3,2)) # Add dummy column for shape
    D = rep("ND", n_in)
    #         Lower Limit          Upper Limit
    i_D =    (y_in[2,] >= 0.98)
    D[i_D] = "PD"
    i_D =                         (y_in[1,] <  0.20) &
                                  (y_in[2,] <  0.98) &
                                  (y_in[3,] <  0.04)
    D[i_D] = "T1"
    i_D =    (y_in[1,] >= 0.20) & (y_in[1,] <  0.50) &
                                  (y_in[3,] <  0.04)
    D[i_D] = "T2"
    i_D =    (y_in[1,] >= 0.50) &
                                  (y_in[3,] <  0.15)
    D[i_D] = "T3"
    i_D = (                       (y_in[1,] <  0.50) &
             (y_in[3,] >= 0.04) & (y_in[3,] <  0.13)) |
          (  (y_in[1,] >= 0.40) & (y_in[1,] <  0.50) &
             (y_in[3,] >= 0.13) & (y_in[3,] <  0.29)) |
          (  (y_in[1,] >= 0.50) &
             (y_in[3,] >= 0.15) & (y_in[3,] <  0.29))
    D[i_D] = "DT"
    i_D =                         (y_in[1,] <  0.23) &
             (y_in[3,] >= 0.13)
    D[i_D] = "D1"
    i_D = (  (y_in[1,] >= 0.23) &
             (y_in[3,] >= 0.29)                    ) |
          (  (y_in[1,] >= 0.23) & (y_in[1,] <  0.40) &
             (y_in[3,] >= 0.13) & (y_in[3,] <  0.29))
    D[i_D] = "D2"
    if(n_in == 1) return(D[1])
    return(D)
  }

  Interpret_DT4 <- function(y_in, n_in) {
    # Almost as in C57.104-2019 Table D.3. I believe there is an error in it.
    # Changed C: (CH4 >= 36 & C2H6 >= 24) to (CH4 >= 36 & C2H6 < 24)
    if(n_in == 1) y_in <- array(rep(unlist(y_in), 2), dim=c(3,2)) # Add dummy column for shape
    D = rep("ND", n_in)
    #         Lower Limit          Upper Limit
    i_D =    (y_in[1,] >= 0.02) & (y_in[1,] <  0.15) &
                                  (y_in[3,] <  0.01)
    D[i_D] = "PD"
    i_D = (  (y_in[2,] >= 0.09) &
             (y_in[3,] >= 0.30) & (y_in[3,] <  0.46)) |
          (  (y_in[2,] >= 0.15) &
             (y_in[3,] >= 0.24) & (y_in[3,] <  0.30)) |
          (                       (y_in[1,] <  0.36) &
             (y_in[3,] >= 0.01) & (y_in[3,] <  0.24)) |
          (  (y_in[1,] >= 0.15) & (y_in[1,] <  0.36) &
             (y_in[3,] <  0.01)) |
          (  (y_in[1,] <  0.02) &
                                  (y_in[3,] <  0.01))
    D[i_D] = "S"
    i_D =                         (y_in[2,] <  0.09) &
             (y_in[3,] >= 0.30)
    D[i_D] = "O"
    i_D = (  (y_in[1,] >= 0.36) &
                                  (y_in[3,] <  0.24)) | # Changed from D.3
          (                       (y_in[2,] <  0.15) &
             (y_in[3,] >= 0.24) & (y_in[3,] <  0.30))
    D[i_D] = "C"
    if(n_in == 1) return(D[1])
    return(D)
  }

  Interpret_DT5 <- function(y_in, n_in) {
    # As in C57.104-2019 Table D.4
    if(n_in == 1) y_in <- array(rep(unlist(y_in), 2), dim=c(3,2)) # Add dummy column for shape
    D = rep("ND", n_in)
    #         Lower Limit          Upper Limit
    i_D =                         (y_in[1,] <  0.01)  &
             (y_in[3,] >= 0.02) & (y_in[3,] <  0.14)
    D[i_D] = "PD"
    i_D = (  (y_in[1,] >= 0.01) & (y_in[1,] <  0.10) &
             (y_in[3,] >= 0.02) & (y_in[3,] <  0.14)) |
          (                       (y_in[1,] <  0.01) &
                                  (y_in[3,] <  0.02)) |
          (                       (y_in[1,] <  0.10) &
             (y_in[3,] >= 0.54))
    D[i_D] = "O"
    i_D =                         (y_in[1,] <  0.10) &
             (y_in[3,] >= 0.14) & (y_in[3,] <  0.54)
    D[i_D] = "S"
    i_D =    (y_in[1,] >= 0.10) & (y_in[1,] <  0.35) &
                                  (y_in[3,] <  0.12)
    D[i_D] = "T2"
    i_D = (  (y_in[1,] >= 0.35) &
                                  (y_in[3,] <  0.12)) |
          (  (y_in[1,] >= 0.50) &
             (y_in[3,] >= 0.12) & (y_in[3,] <  0.14)) |
          (  (y_in[1,] >= 0.70) &
             (y_in[3,] >= 0.14)) |
          (  (y_in[1,] >= 0.35) &
             (y_in[3,] >= 0.30))
    D[i_D] = "T3"
    i_D = (  (y_in[1,] >= 0.10) & (y_in[1,] <  0.50) &
             (y_in[3,] >= 0.12) & (y_in[3,] <  0.14)) |
          (  (y_in[1,] >= 0.10) & (y_in[1,] <  0.70) &
             (y_in[3,] >= 0.14) & (y_in[3,] <  0.30))
    D[i_D] = "C"
    if(n_in == 1) return(D[1])
    return(D)
  }

  if(DT_n == "1") return(Interpret_DT1(y_in, n_in))
  if(DT_n == "4") return(Interpret_DT4(y_in, n_in))
  if(DT_n == "5") return(Interpret_DT5(y_in, n_in))
}
