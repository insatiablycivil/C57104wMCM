#' Interpret_IECn: Applies IEC Table 1 or Table 2 Diagnostic Logic
#'
#' This functions interprets the inputted gas sample ratios to determine the appropriate diagnosis for each in accordance to either IEC 60599:2022 Table 1 or Table 2. For Table 1, bounds shown in Figure B.1 are used.
#'
#' @param IEC_n String: Value corresponding to either Table 1 or Table 2. Expected value: 1 or 2.
#' @param y_in Numeric Array: Values representing generated gas ratios, assumed order: C2H2/C2H4, CH4/H2, C2H4/C2H6, assumed shape: 3 x Samples.
#' @param n_in Numeric: Number of samples to be diagnosed.
#' @return Strings: List of Diagnoses, one for each inputted sample.
#' @examples
#' Interpret_IECn(IEC_n = "1",
#'                y_in = array(c(0.090, 0.002, 0.401, 0.120, 0.001, 0.310), dim = c(3,2)),
#'                n_in = 2)
#' @export
Interpret_IECn <- function(IEC_n, y_in, n_in) {
  # Assumed y_n order: C2H2/C2H4, CH4/H2, C2H4/C2H6
  Interpret_IEC1 <- function(y_in, n_in) {
    # Bounds based on IEC 60599:2022 Table 1
    # Table uses "NS" defined as "Non-significant whatever the value".
    # This phrasing is bad: to me, that means any value.
    # Figure B.1 interprets it as 0.01, but then for C2H2/C2H4, it has "> 1 but NS"..?
    # Using definitions based on Figure B.1, despite my dislike for it.
    if(n_in == 1) y_in <- array(rep(unlist(y_in), 2), dim=c(3,2)) # Add dummy column for shape
    D = rep("ND", n_in)
    #         Lower Limit          Upper Limit
    i_D =                         (y_in[1,] <  0.01) & # See note on "NS"
                                  (y_in[2,] <  0.1) &
                                  (y_in[3,] <  0.2)
    D[i_D] = "PD"
    i_D =                         (y_in[1,] <  0.01) & # See note on "NS"
              (y_in[2,] >  1.0) &                      # See note on "> 1 but NS"
                                  (y_in[3,] <  1.0)
    D[i_D] = "T1"
    i_D =                         (y_in[1,] <  0.1) &
              (y_in[2,] >  1.0) &
              (y_in[3,] >= 1.0) & (y_in[3,] <= 4.0)
    D[i_D] = "T2"
    i_D =                         (y_in[1,] <  0.2) &
              (y_in[2,] >  1.0) &
                                  (y_in[3,] >  4.0)
    D[i_D] = "T3"
    i_D =     (y_in[1,] >= 0.6) & (y_in[1,] <= 2.5) &
              (y_in[2,] >= 0.1) & (y_in[2,] <= 1.0) &
              (y_in[3,] >  2.0)
    D[i_D] = "D2"
    i_D =     (y_in[1,] >  1.0) &
              (y_in[2,] >= 0.1) & (y_in[2,] <= 0.5) &
              (y_in[3,] >  1.0)
    D[i_D] = "D1"
    if(n_in == 1) return(D[1])
    return(D)
  }

  Interpret_IEC2 <- function(y_in, n_in) {
    # Bounds based on IEC 60599:2022 Table 2
    if(n_in == 1) y_in <- array(rep(unlist(y_in), 2), dim=c(3,2)) # Add dummy column for shape
    D = rep("ND", n_in)
    #         Lower Limit          Upper Limit
    i_PD =                        (y_in[2,] <  0.2)
    D[i_PD] = "PD"
    i_D =                         (y_in[1,] <  0.2)
    D[i_D] = "T"
    if (sum(i_D & i_PD)>0) warning("Warning: IEC Table 2 could be both PD and T, T chosen!")
    i_D =    (y_in[1,] >  0.2)
    D[i_D] = "D"
    if (sum(i_D & i_PD)>0) warning("Warning: IEC Table 2 could be both PD and D, D chosen!")
    if(n_in == 1) return(D[1])
    return(D)
  }

  if(IEC_n == "1") return(Interpret_IEC1(y_in, n_in))
  if(IEC_n == "2") return(Interpret_IEC2(y_in, n_in))
}
