#' Preprocess: Validate inputs and preprocess them into expected data formats / values
#'
#' This function is intended to capture errors. No further validation is then done beyond this function.
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
#' @param silent Boolean: Defaulted to FALSE.
#' @param verbose Boolean: Defaulted to FALSE.
#' @param validate Boolean: Defaulted to FALSE.
#' @return List containing:
#' \itemize{
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
#'  \item $\code{sd_data}: Data.frame. Same as \code{gas_data} but with the values representing standard deviation and not gas values.
#' }
#' @examples
#' Preprocess(gas_data = data.frame(date = as.Date(c("2024-06-20","2024-06-21","2024-06-22","2024-06-23","2024-06-24")), H2 = c(8,4,2,1,0.5), CH4 = c(10,20,40,80,160)),
#'            limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
#'            accuracy = setNames(c(0.15, 0.50), c("H2", "CH4")),
#'            min_accuracy = setNames(c(1, 5), c("H2", "CH4")),
#'            distribution_shape = "N",
#'            rho = 0.5,
#'            diagnosis = 0,
#'            n = 1000,
#'            flex = FALSE, silent = FALSE, verbose = FALSE, validate = FALSE)
#' Preprocess(gas_data = data.frame(date=as.Date(c("2024-06-20","2024-06-21","2024-06-22","2024-06-23","2024-06-24")), H2=c(80,80,80,80,82),CH4=c(90,95,100,105,107),C2H6=c(80,100,80,80,82),C2H4=c(30,30,31,29,31),C2H2=c(1,1,0,1,3),CO=c(800,850,860,850,852),CO2=c(900,950,940,940,942)),
#'            limit_data = data.frame(H2=c(80,200,40,20),CH4=c(90,150,30,10),C2H6=c(90,175,25,9),C2H4=c(50,100,20,7),C2H2=c(1,2,0.5,0.5),CO=c(900,1100,20,100),CO2=c(9000,12500,250,1000)),
#'            accuracy = 0.15,
#'            min_accuracy = 1,
#'            distribution_shape = "T",
#'            rho = matrix(rep(1,49), ncol=7),
#'            diagnosis = "DT",
#'            n = 1000,
#'            flex = TRUE, silent = FALSE, verbose = TRUE, validate = TRUE)
#' @export
Preprocess <- function(gas_data, limit_data, accuracy, min_accuracy,
                       distribution_shape, rho, diagnosis, n,
                       flex = FALSE, silent = FALSE, verbose = FALSE, validate = FALSE) {
  # Preprocess will validate inputs and convert to expected format
  gas_names = c("H2", "CH4", "C2H6", "C2H4", "C2H2", "CO", "CO2")

  Check_Wordage <- function(value, word) {
    if(is.null(value)) stop(paste0("Error: \'", word,"\' must be a boolean equivalent"))
    if(is.na(as.logical(value))) stop(paste0("Error: \'", word,"\' must be a boolean equivalent"))
    return(as.logical(value))
  }

  Check_Samples <- function(n) {
    if(is.null(n)) stop("Error: \'n\' must be a numeric equivalent")
    if(is.na(as.numeric(n))) stop("Error: \'n\' must be a numeric equivalent")
    if(n<=0) stop("Error: \'n'\ must be above 0")
    if(n==1) stop("Error: \'n'\ cannot be 1")
    if((n>1)&(n<10)) warning("Warning: \'n\' is very low, results will be unreliable and unstable!")
    if(n<1){ # Treat as error tolerance
      k = 1.96 # Confidence Probability of 95%
      P = 0.5  # Maximum Uncertainty when at 50% Probability of success
      # Rearrangement of Wald Interval, 'n' is interpreted as W
      return(ceiling((P*(1-P))/((n)/(2*k))^2))
    } else return(n) # Treat as number of samples
  }

  Check_Accuracy <- function(U, V) {
    if(is.null(U)) return(0)
    if((length(U)==1) && is.na(U)) return(0)
    if(V=="accuracy") {
      if(any(U<0) | any(U>=1)) stop("Error: \'accuracy\' must be between 0-1, exclusive (use NA to exclude)")
    } else if(any(U<0)) stop("Error: \'min_accuracy' must be positive (use NA to exclude)")
    if((length(U)>1)) {
      if(any(is.na(U))) stop(paste0("Error: If using a \'", V, "\' per gas, cannot contain NAs"))
      if(is.null(names(U))) stop(paste0("Error: If using a \'", V, "\' per gas, must used named variable"))
    }
    return(U)
  }

  Check_Shape <- function(distribution_shape) {
    if(distribution_shape %in% c("N","n","Norm","norm","Normal","normal")) return(1)
    if(distribution_shape %in% c("T","t","Tri","tri","Triangle","triangle")) {
      library(triangle)
      return(2)
    }
    if(distribution_shape %in% c("U","u","Uni","uni","Uniform","uniform")) return(3)
    stop("Error: \'distribution_shape\' must be either: <uniform|triangle|normal>")
  }

  Check_Diagnosis <- function(diagnosis, gas_data) {
    if(is.numeric(diagnosis) && (diagnosis==0)) return(diagnosis)
    if(!all(c("H2", "CH4", "C2H6", "C2H4", "C2H2") %in% colnames(gas_data))) {
      stop("Error: \'gas_data\' must have the following gases for diagnosis: \n<H2,CH4,C2H6,C2H4,C2H2,CO,CO2>")
    }
    if(is.numeric(diagnosis)) {
      if(diagnosis==1) return(diagnosis)
    }
    if(diagnosis == "DT") return(2)
    if(diagnosis == "RR") return(3)
    if(diagnosis == "IEC") return(4)
    stop("Error: \'diagnosis'\ must be either: <0|1> or <DT|RR|IEC>")
  }

  Check_Rho <- function(rho, flex) {
    if(is.null(rho)) return(0)
    if(length(rho)==1) { # Assumed rho
      if(is.na(rho)) return(0)
      if((rho<0)|(rho>1)) stop("Error: Assumed \'rho\' must be between 0-1, exclusive (use NA to exlude)")
      return(rho)
    } # Assumed cor_matrix
    if(is.null(dim(rho))) stop("Error: \'rho\' assumed correlation matrix and must be in a matrix-like form, (having dim())")
    if(any(is.na(rho))) stop("Error: \'rho\' assumed correlation matrix and cannot contain NAs (use single NA to exclude)")
    library(corpcor)
    if(!is.positive.definite(rho)) {
      if(!flex) stop("Error: \'rho\' assumed correlation matrix is not positive definite, use flex == TRUE")
      warning("Warning: \'rho\' assumed correlation matrix and changed into positive definite")
      rho <- make.positive.definite(rho)
      print(rho)
    }
    return(rho)
  }

  Check_Gas <- function(gas_data) {
    is.date <- function(x) inherits(x, 'Date')
    if(!inherits(gas_data,"data.frame")){
      stop("Error: gas_data must be a data.frame and each row being new sample date (YYYY-MM-DD).")
    }
    if(!any(colnames(gas_data) %in% gas_names)) {
      stop("Error: gas_data must have at least one of these column names: \n<H2,CH4,C2H6,C2H4,C2H2,CO,CO2>")
    }
    if(!all(sapply(gas_data$date, is.date))){
      stop("Error: gas_data must have a column 'date' with dates (YYYY-MM-DD).")
    }
    # Want them ordered with newest at the top, but listed as days since oldest.
    if(nrow(gas_data)==1) {
      gas_data$date <- 0
      return(gas_data)
    }
    if(nrow(gas_data)>6){
      warning("Warning: gas_data has more rows (>6) than expected for Table 4! All will be used.")
    }
    mask <- rev(order(gas_data$date))
    gas_data <- gas_data[mask,]
    if(all(diff(gas_data$date)<0)) {
      gas_data$date <- as.numeric(gas_data$date - gas_data$date[nrow(gas_data)])
    } else stop("Error: gas_data cannot have tests taken on same day")
    return(gas_data)
  }

  Check_Limit <- function(limit_data) {
    # Not bothering to check it's numeric
    if(!inherits(limit_data,"data.frame")){
      stop("Error: limit_data must be a data.frame and each row being Table Limits [PPM-equiv].")
    }
    if(!all(colnames(limit_data) %in% gas_names)) {
      stop("Error: limit_data must have columns for all gases included in gas_data.")
    }
    if(nrow(limit_data) != 4){
      stop("Error: limit_data must have 4 rows for Table Limits 1-4 [PPM-equiv], inclusive.")
    }
    limit_data[4,] <- limit_data[4,]/365.25 # Convert to per day rate
    return(limit_data)
  }

  ### Validation
  silent = Check_Wordage(silent, "silent")
  verbose  = Check_Wordage(verbose, "verbose")
  validate  = Check_Wordage(validate, "validate")
  n = Check_Samples(n)
  accuracy = Check_Accuracy(accuracy, "accuracy")
  min_accuracy = Check_Accuracy(min_accuracy, "minimum accuracy")
  if(length(accuracy) != length(min_accuracy)) stop("Error: accuracy and min_accuracy must be of same length")
  if((length(accuracy) > 1) &&
     any(is.na(match(names(accuracy), names(min_accuracy))))) stop("Error: accuracy and min_accuracy must have same labels")
  distribution_shape = Check_Shape(distribution_shape)
  diagnosis = Check_Diagnosis(diagnosis, gas_data)
  gas_data = Check_Gas(gas_data)
  limit_data = Check_Limit(limit_data)
  rho = Check_Rho(rho, flex) # Too onerous to robustly validate this: assumed order correct

  # Reordering gases and limits
  gas_data <- gas_data[,c('date',intersect(gas_names,colnames(gas_data)))]
  col_names <- colnames(gas_data)[-1]
  if(ncol(limit_data) != length(col_names)) {
    col_names <- intersect(col_names, colnames(limit_data))
    warning("Warning: Potential column-mismatch between \'gas_data\' and \'gas_limits\'.")
    limit_data <- limit_data[,col_names]
    gas_data <- gas_data[,c('date', col_names)]
  }

  if((length(accuracy) > 1)) {
    if(length(accuracy) != length(col_names)) warning("Warning: Potential column-mismatch between \'gas_data\' and \'accuracy\' and/or \'min_accuracy\'")
    mask <- names(accuracy) %in% col_names
    if(sum(mask) != length(col_names)) stop("Error: \'accuracy\' is missing values for some expected gases")
    accuracy <- accuracy[mask]
    mask <- names(min_accuracy) %in% col_names
    min_accuracy <- min_accuracy[mask]
    min_accuracy <- min_accuracy[match(names(accuracy), names(min_accuracy))]
  }

  if((length(col_names) == 1) &&
     (any(!(is.na(rho))) && any(rho != 0))) {
    warning("Warning: Only single gas present, \'rho\' will be ignored.")
    rho <- 0
    cor_mat <- NA
  } else {
    if(length(rho)>1) {
      warning("Warning: \'rho\' assumed correlation matrix: ensure columns ordered: <H2,CH4,C2H6,C2H4,C2H2,CO,CO2>")
      if(ncol(rho) != length(col_names)) stop("Error: \'rho\' assumed correlation matrix but does not match number of gases")
      cor_mat <- rho
      rho <- NA
    } else {
      cor_mat <- NA
    }
  }

  # Separate Dates and Gases
  sd_data <- temp_acc <- temp_min_acc <- gas_data # Gas values will be replaced with SD values.
  temp_acc[,col_names] <- rep(rep(accuracy,length(col_names)/length(accuracy)),each= nrow(gas_data)) # Reshaping to match sd_data
  temp_min_acc[,col_names] <- rep(rep(min_accuracy,length(col_names)/length(min_accuracy)),each= nrow(gas_data)) # Reshaping to match sd_data
  sd_data[,col_names] <- switch(distribution_shape,
                                pmax(temp_acc[,col_names]*sd_data[,col_names], temp_min_acc[,col_names])/1.96,        # Norm
                                sqrt((pmax(temp_acc[,col_names]*sd_data[,col_names], temp_min_acc[,col_names])^2)/6), # Tri
                                sqrt((pmax(temp_acc[,col_names]*sd_data[,col_names], temp_min_acc[,col_names])^2)/3), # Uni
                                stop("Error: Unknown \'distribution_shape\'!"))

  output <- list(gas_data = gas_data,
                 limit_data = limit_data,
                 accuracy = accuracy,
                 min_accuracy = min_accuracy,
                 distribution_shape = distribution_shape,
                 rho = rho,
                 diagnosis = diagnosis,
                 n = n,
                 col_names = col_names,
                 cor_mat = cor_mat,
                 sd_data = sd_data)

  return(output)
}
