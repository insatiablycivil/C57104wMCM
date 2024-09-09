#' C57104wMCM: Simplified C57.104 DGA Screening with propagated uncertainty via a Monte Carlo Method
#'
#' This function applies a simplified version of the C57.104 DGA Screening methodology with propagated uncertainty. This is achieved by generating random samples according to the statistical description. These outcomes for these random samples are used as proxies for estimating the likelihood of different outcomes. There is also an integrated diagnostic step available using a combination of Duval Triangles 1|4|5, Rogers Ratio Method and IEC 60599:2022 Tables 1|2.
#'
#' @param gas_data Data.frame: First column containing date of sample. Cannot have two samples on same day. Rest of columns containing sample value in PPM.
#' @param limit_data Data.frame: Columns containing the IEEE DGA Table limits for each gas. Each Table is on its own row; four rows expected.
#' @param accuracy Named Numeric: Values to be multiplied by gas value to give uncertainty. Can either be same value for all gases, or individual values for each gas. Expected range: 0-1.
#' @param min_accuracy Named Numeric: Absolute Values (PPM) used as minimum uncertainty. Can either be same value for all gases, or individual values for each gas.
#' @param distribution_shape Numeric: Distribution shape representing uncertainty. "N" = Normal, "T" = Triangular, "U" = Uniform. If Normal, Accuracy is interpreted as 95% confidence interval, else as 100%.
#' @param rho Numeric: Value representing rho for correlations. Can either be same value for all gases, individual values for each gas, or a matrix intended to represent the \code{cor_mat}.
#' @param diagnosis Mixed: Used to specify which diagnoses to conduct. 0 = None, 1 = All, "DT" = Duval Triangles 1|4|5 only, "RR" = Rogers Ratio Only, "IEC" = IEC Tables 1|2 Only.
#' @param n Numeric: Number of samples to generate for the Monte Carlo for each gas. If a value <1, it is assumed as an estimate for the error tolerance, then Wald Interval used to estimate n.
#' @param flex Boolean: Defaulted to FALSE. If TRUE, and \code{rho} is assumed a correlation matrix equivalent, it will attempt to force it to become a positive definite if needed.
#' @param silent Boolean: Defaulted to FALSE. If TRUE: heavily restricts outputs to console.
#' @param verbose Boolean: Defaulted to FALSE. If FALSE: restricts outputs to console.
#' @param validate Boolean: Defaulted to FALSE. If TRUE: provides various validatory plots.
#' @return Nested List containing:
#' \itemize{
#'  \item $\code{inputs}:
#'  \itemize{
#'  \item $\code{gas_data}: Data.frame. First column as Date converted to days since oldest sample. Rest as gas values in PPM.
#'  \item $\code{limit_data}: Data.frame. Unchanged from input.
#'  \item $\code{accuracy}: Numeric. Unchanged from input.
#'  \item $\code{min_accuracy}: Numeric. Unchanged from input.
#'  \item $\code{distribution_shape}: Numeric. 1 = Normal, 2 = Triangular, 3 = Uniform.
#'  \item $\code{rho}: Numeric or NA. NA if \code{cor_mat} is being used.
#'  \item $\code{diagnosis}: Numeric. 0 = None, 1 = All, 2 = Duval Triangles, 3 = Rogers Ratios, 4 = IEC Tables.
#'  \item $\code{n}: Numeric. Number of samples used for Monte Carlo Simulation.
#'  \item $\code{col_names}: Names of found gas names.
#'  \item $\code{cor_mat}: Matrix or NA. NA if \code{rho} is being used.
#'  \item $\code{sd_data}: Data.frame. Same as \code{gas_data} but with the values representing Std Dev and not gas values.
#'  }
#'  \item $\code{sampled}: Optional: (\code{validate} == TRUE): Numeric Array. Values representing generated samples for the Monte Carlo based on the statistical description provided. Shaped: Days x Gases x Samples.
#'  \item $\code{screened}:
#'  \itemize{
#'  \item $\code{protocol}: Numeric. 1|2|3 depending on if 1, 2, or more samples.
#'  \item $\code{T<n>_default}: Output for Tables 1|2|3|4 using default gas values.
#'  \item $\code{T<n>}: Probability of output for Tables 1|2|3|4 using generated samples.
#'  \item $\code{T<n>_bounds}: Numeric. Associated estimate of confidence interval.
#'  \item $\code{gas_L<n>_default}: Numeric. Output for per-gas DGA Status Levels 1|2|3.
#'  \item $\code{gas_L<n>}: Numeric. Probability of output for per-gas DGA Status Levels 1|2|3.
#'  \item $\code{gas_L<n>_bounds}: Numeric. Associated confidence interval.
#'  \item $\code{L<n>_default}: Numeric. Output for combined DGA Status Levels 1|2|3.
#'  \item $\code{L<n>}: Numeric. Probability of output for combined DGA Status Levels 1|2|3.
#'  \item $\code{L<n>_bounds}: Numeric. Mask indicated which generated samples were combined DGA Status Levels 1|2|3.
#'  \item $\code{L<n>_mask}: Optional: (\code{validate} == TRUE): Logical Array. Associated estimate of confidence interval.
#'  \item $\code{T12_vals_default}: Optional: (\code{validate} == TRUE): Numeric. Metric values using default gas values used for Tables 1 and 2.
#'  \item $\code{T3_vals_default}: Optional: (\code{validate} == TRUE): Numeric. Metric values using default gas values used for Table 3.
#'  \item $\code{T4_vals_default}: Optional: (\code{validate} == TRUE): Numeric. Metric values using default gas values used for Tables 4.
#'  \item $\code{T12_vals}: Optional: (\code{validate} == TRUE): Numeric. Metric values using generated samples used for Tables 1 and 2.
#'  \item $\code{T3_vals}: Optional: (\code{validate} == TRUE): Numeric. Metric values using generated samples used for Table 3.
#'  \item $\code{T4_vals}: Optional: (\code{validate} == TRUE): Numeric. Metric values using generated samples used for Tables 4.
#'  \item $\code{GL<n>_mask}: Optional: (\code{validate} == TRUE): Logical Array. Mask indicated which generated samples were per-gas DGA Status Levels 1|2|3.
#'  }
#'  \item $\code{diagnoses}: Optional: (\code{diagnosis} != 0):
#'  \itemize{
#'  \item $\code{DT<n>_default_data}: Optional: (\code{diagnosis} == 1|2): Numeric. Ratios relevant to DT 1|4|5 using default gas values.
#'  \item $\code{DT<n>_0_data}: Optional: (\code{diagnosis} == 1|2): Numeric Array. Ratios relevant to DT 1|4|5 using generated samples.
#'  \item $\code{DT<n>_L<n>}: Optional: (\code{diagnosis} == 1|2): Numeric. Proportion of samples, screened at DGA Status 1|2|3, assigned given diagnosis in accordance to DT 1|4|5.
#'  \item $\code{DT<n>_Ln1}: Optional: (\code{diagnosis} == 1|2): Numeric. Proportion of samples, screened at DGA Status 2 or 3, assigned given diagnosis in accordance to DT 1|4|5.
#'  \item $\code{gas_L<n>_default}: Optional: (\code{diagnosis} == 1|3|4): Numeric. Output for per-gas DGA Status Levels 1|2|3.
#'  \item $\code{Abs_default_default}: Optional: (\code{diagnosis} == 1|3|4): Numeric. Ratios relevant to RR and/or IEC 1|2 using default gas values.
#'  \item $\code{Abs_0_data}: Optional: (\code{diagnosis} == 1|3|4): Numeric Array. Ratios relevant to RR and/or IEC 1|2 using generated samples.
#'  \item $\code{RR_L<n>}: Optional: (\code{diagnosis} == 1|3): Numeric. Proportion of samples, screened at DGA Status 1|2|3, assigned given diagnosis in accordance to RR.
#'  \item $\code{RR_Ln1}: Optional: (\code{diagnosis} == 1|3): Numeric. Proportion of samples, screened at DGA Status 2 or 3, assigned given diagnosis in accordance to RR.
#'  \item $\code{IEC<n>_L<n>}: Optional: (\code{diagnosis} == 1|4): Numeric. Proportion of samples, screened at DGA Status 1|2|3, assigned given diagnosis in accordance to IEC 1|2.
#'  \item $\code{IEC<n>_Ln1}: Optional: (\code{diagnosis} == 1|4): Numeric. Proportion of samples, screened at DGA Status 2 or 3, assigned given diagnosis in accordance to IEC 1|2.
#'  }
#' }
#' @examples
#' # Example 1: Two gases, two samples, individually assigned accuracies, no correlation, and no diagnostics
#' # gas_data: Dates and gas values (PPM). With only two samples, Table 4 will not be used.
#' gas_data_in <- data.frame(date = as.Date(c("2024-06-20","2024-06-24")),
#'                           H2 = c(8,4), CH4 = c(80,160))
#' # limit_data: IEEE Tables 1-4 Limits for each gas (PPM and PPM/Year).
#' limit_data_in <- data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03))
#' # accuracy: Multiplied by gas value to give uncertainty.
#' # If "Normal", this represents 95% probability level, else 100%.
#' # Can also give a single value to be used for all gases
#' accuracy_in <- setNames(c(0.15, 0.50), c("H2", "CH4"))
#' # min_accuracy: Treated as PPM and used if the above would give an otherwise smaller value.
#' # Either can be given as NA to force the use of the other. This too can be given as single value.
#' min_accuracy_in = setNames(c(1, 5), c("H2", "CH4"))
#' C57104wMCM(gas_data = gas_data_in,
#'            limit_data = limit_data_in,
#'            accuracy = accuracy_in,
#'            min_accuracy = min_accuracy_in,
#'            distribution_shape = "N",
#'            rho = 0.5,
#'            diagnosis = 0,
#'            n = 1000,
#'            flex = FALSE, silent = FALSE, verbose = FALSE, validate = FALSE)
#'
#' # Example 2: Seven gases, five samples, shared accuracies, specified correlation matrix, and Duval Triangle diagnostics
#' gas_data_in <- data.frame(date=as.Date(c("2024-06-20","2024-06-21","2024-06-22","2024-06-23","2024-06-24")),
#'                           H2=c(80,80,80,80,82), CH4=c(90,95,100,105,107), C2H6=c(80,100,80,80,82),
#'                           C2H4=c(30,30,31,29,31), C2H2=c(1,1,0,1,3), CO=c(800,850,860,850,852),
#'                           CO2=c(900,950,940,940,942))
#' limit_data_in <- data.frame(H2=c(80,200,40,20),CH4=c(90,150,30,10),C2H6=c(90,175,25,9),
#'                             C2H4=c(50,100,20,7),C2H2=c(1,2,0.5,0.5),CO=c(900,1100,20,100),
#'                             CO2=c(9000,12500,250,1000))
#' # rho: If given as a matrix, it is interpreted as a correlation matrix, else as rho.
#' # In this example, the correlation matrix is not positive definite.
#' # However, the flex argument allows for it to be adjusted to make it so.
#' rho_in <- matrix(rep(1,49), ncol=7)
#' C57104wMCM(gas_data = gas_data_in,
#'            limit_data = limit_data_in,
#'            accuracy = 0.15,
#'            min_accuracy = 1,
#'            distribution_shape = "T",
#'            rho = rho_in,
#'            diagnosis = "DT",
#'            n = 500,
#'            flex = TRUE, silent = FALSE, verbose = TRUE, validate = TRUE)
#' @export
C57104wMCM <- function(gas_data, limit_data, accuracy, min_accuracy,
                       distribution_shape = "n", rho = 0, diagnosis = 0, n = 1000,
                       flex = FALSE, silent = FALSE, verbose = FALSE, validate = FALSE) {
  # Preprocess() will validate inputs and convert to expected format
  # Sample() will generate relevant samples from stated distributions
  # Screening() will obtain IEEE Screening outputs for generated samples
  # Diagnosis() will obtain diagnoses
  output <- list()

  preprocessed <- C57104wMCM::Preprocess(gas_data, limit_data, accuracy, min_accuracy,
                                         distribution_shape, rho, diagnosis, n,
                                         flex, silent, verbose, validate)
  if(is.character(preprocessed)) return(preprocessed)
  output$inputs <- preprocessed

  # output is [days, gases, samples]
  y_N <- C57104wMCM::Sample(preprocessed$gas_data, preprocessed$limit_data,
                            preprocessed$accuracy, preprocessed$min_accuracy,
                            preprocessed$distribution_shape, preprocessed$rho,
                            preprocessed$cor_mat, preprocessed$n, preprocessed$col_names,
                            preprocessed$sd_data)
  if(validate == TRUE) {
    output$sampled <- y_N
    C57104wMCM::Sample_Plot(y_N, preprocessed)
  }

  screened <- C57104wMCM::Screening(y_N, preprocessed$limit_data, preprocessed$gas_data, silent, verbose)
  if(validate == TRUE) {
    output$screening <- screened
    C57104wMCM::Screening_Plot(screened, preprocessed)
  } else output$screening <- screened[1:31] # Excluding masks and raw metric values

  if(preprocessed$diagnosis == 0) return(output)
  diagnosed <- C57104wMCM::Diagnosis(y_N[1,,], preprocessed$gas_data, preprocessed$diagnosis, preprocessed$col_names,
                                     screened$L1_mask, screened$L2_mask, screened$L3_mask, silent, verbose)
  output$diagnoses = diagnosed
  if(validate == TRUE) C57104wMCM::Diagnosis_Plot(diagnosed, screened, preprocessed)
  return(output)
}
