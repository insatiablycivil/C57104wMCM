#' Diagnosis: Apply Diagnostics in accordance to either Duval Triangles, Rogers Ratios, or IEC Tables, or all.
#'
#' This function applies diagnostics in accordance to either Duval Triangles, Rogers Ratios, or IEC Tables, or all. The results are separated based on the IEEE Screening DGA Status Levels.
#'
#' @param y_1 Numeric Array: Values representing generated samples for the Monte Carlo based on the statistical description provided. Shaped: Gases x Samples.
#' @param gas_data Data.frame: First column containing days since oldest sample. Rest of columns containing sample value in PPM.
#' @param diagnosis Data.frame: Columns containing the IEEE DGA Table limits for each gas. Each Table is on its own row; four rows expected. Caution: Table 4 is expressed as PPM/Day!
#' @param col_names Strings: Gas Names. Order must match other named variables throughout. Expected a subset of <H2|CH4|C2H6|C2H4|C2H2|CO|CO2>.
#' @param i_L1 Logical Array. Mask indicated which generated samples were per-gas DGA Status Level 1. Length equal to n samples.
#' @param i_L2 Logical Array. Mask indicated which generated samples were per-gas DGA Status Level 2. Length equal to n samples.
#' @param i_L3 Logical Array. Mask indicated which generated samples were per-gas DGA Status Level 3. Length equal to n samples.
#' @param silent Boolean: Defaulted to FALSE. If TRUE heavily restricts outputs to console.
#' @param verbose Boolean: Defaulted to FALSE. If FALSE restricts outputs to console.
#' @return List containing a subset of the following:
#' \itemize{
#' \item $\code{DT<n>_default_data}: Optional: (\code{diagnosis} == 1|2): Numeric. Ratios relevant to DT 1|4|5 using default gas values.
#' \item $\code{DT<n>_0_data}: Optional: (\code{diagnosis} == 1|2): Numeric Array. Ratios relevant to DT 1|4|5 using generated samples.
#' \item $\code{DT<n>_L<n>}: Optional: (\code{diagnosis} == 1|2): Numeric. Proportion of samples, screened at DGA Status 1|2|3, assigned given diagnosis in accordance to DT 1|4|5.
#' \item $\code{DT<n>_Ln1}: Optional: (\code{diagnosis} == 1|2): Numeric. Proportion of samples, screened at DGA Status 2 or 3, assigned given diagnosis in accordance to DT 1|4|5.
#' \item $\code{gas_L<n>_default}: Optional: (\code{diagnosis} == 1|3|4): Numeric. Output for per-gas DGA Status Levels 1|2|3.
#' \item $\code{Abs_default_default}: Optional: (\code{diagnosis} == 1|3|4): Numeric. Ratios relevant to RR and/or IEC 1|2 using default gas values.
#' \item $\code{Abs_0_data}: Optional: (\code{diagnosis} == 1|3|4): Numeric Array. Ratios relevant to RR and/or IEC 1|2 using generated samples.
#' \item $\code{RR_L<n>}: Optional: (\code{diagnosis} == 1|3): Numeric. Proportion of samples, screened at DGA Status 1|2|3, assigned given diagnosis in accordance to RR.
#' \item $\code{RR_Ln1}: Optional: (\code{diagnosis} == 1|3): Numeric. Proportion of samples, screened at DGA Status 2 or 3, assigned given diagnosis in accordance to RR.
#' \item $\code{IEC<n>_L<n>}: Optional: (\code{diagnosis} == 1|4): Numeric. Proportion of samples, screened at DGA Status 1|2|3, assigned given diagnosis in accordance to IEC 1|2.
#' \item $\code{IEC<n>_Ln1}: Optional: (\code{diagnosis} == 1|4): Numeric. Proportion of samples, screened at DGA Status 2 or 3, assigned given diagnosis in accordance to IEC 1|2.
#' }
#' @examples
#' # y_1: Numeric array of generated samples, shaped: Gases, Samples
#' y_1_in <- array(c(82.39, 107.51,  82.39,  31.15,   3.03, 856.02, 946.45,  77.95, 101.71,  77.95,  29.47,   2.67, 809.89, 895.44,
#'                   78.98, 103.06,  78.98,  29.86,   2.75, 820.60, 907.28,  84.62, 110.42,  84.62,  31.99,   3.21, 879.26, 972.14,
#'                   74.19,  96.81,  74.19,  28.05,   2.37, 770.88, 852.31,  85.56, 111.64,  85.56,  32.35,   3.29, 888.97, 982.88,
#'                   84.71, 110.54,  84.71,  32.03,   3.22, 880.20, 973.18,  74.09,  96.68,  74.09,  28.01,   2.36, 769.80, 851.12,
#'                   83.24, 108.62,  83.24,  31.47,   3.10, 864.91, 956.27,  76.22,  99.45,  76.22,  28.81,   2.53, 791.90, 875.55), dim = c(7,10))
#' gas_data_in <- data.frame(date=c(4,3,2,1,0),
#'                           H2=c(80,80,80,80,82), CH4=c(90,95,100,105,107), C2H6=c(80,100,80,80,82),
#'                           C2H4=c(30,30,31,29,31), C2H2=c(1,1,0,1,3), CO=c(800,850,860,850,852),
#'                           CO2=c(900,950,940,940,942))
#' # temp_rndm_screen is just for generating screening mask as an example.
#' temp_rndm_screen <- sample(1:3, 10, replace=TRUE)
#' i_L1_in <- temp_rndm_screen == 1
#' i_L2_in <- temp_rndm_screen == 2
#' i_L3_in <- temp_rndm_screen == 3
#'
#' Diagnosis(y_1 = y_1_in,
#'           gas_data = gas_data_in,
#'           diagnosis = 1,
#'           col_names = c("H2", "CH4", "C2H6", "C2H4", "C2H2", "CO", "CO2"),
#'           i_L1 = i_L1_in,
#'           i_L2 = i_L2_in,
#'           i_L3 = i_L3_in,
#'           silent = FALSE, verbose = FALSE)
#' @export
Diagnosis <- function(y_1, gas_data, diagnosis, col_names, i_L1, i_L2, i_L3, silent = FALSE, verbose = FALSE) {
  # Diagnosis: <1|2|3|4> == <ALL|DT|RR|IEC>
  # diagnosis_logic_<DTn|IECn|RR> contain functions for diagnosing limits
  # Bunch of functions for printing outputs
  output <- list()

  if(silent == FALSE) {
    print_message <- function(label, val_1, val_2, val_3) {
      label <- sprintf("%5s", label)
      val_1 <- sprintf("%06.2f", val_1)
      val_2 <- sprintf("%06.2f", val_2)
      val_3 <- sprintf("%06.2f", val_3)
      cat(label, ":", val_1, "(", val_2, "/", val_3, ")\n")
    }

    print_rel_ratios <- function(y_in, y_0_in, mask_in, i_L1_in, i_L2_in, i_L3_in, skip) {
      print_template <- function(var_name, i_in, y_0_in) {
        if(sum(i_in)==0) return(cat(paste0("\nNo examples of ", var_name, "!\n")))
        cat(paste0("\nRatios for ", var_name, ":\n"))
        for(idx in c(1:3)) {
          a <- paste0(col_names[mask_in][[idx]],"%")
          b <- sum(y_0_in[idx,i_in])
          N <- length(y_0_in[idx,i_in])
          print_message(a, (b/N)*100, b, N)
        }
      }
      N <- sum(y_in)
      cat("\nDefault Ratios:\n")
      for(idx in c(1:3)) {
        a <- paste0(col_names[mask_in][[idx]],"%")
        b <- y_in[[idx]]
        print_message(a, (b/N)*100, b, N)
      }
      N <- length(i_L1_in)
      if(skip == 0) cat("\nRatios for ALL:\n")
      if(skip != 0) cat(paste0("\nRatios for ALL (All == DGA Status ", skip, "):\n"))
      for(idx in c(1:3)) {
        a <- paste0(col_names[mask_in][[idx]],"%")
        b <- sum(y_0_in[idx,])
        print_message(a, (b/N)*100, b, N)
      }
      if((verbose == TRUE) & (skip == 0)) {
        print_template("SL1", i_L1_in, y_0_in)
        print_template("SL2", i_L2_in, y_0_in)
        print_template("SL3", i_L3_in, y_0_in)
      }
    }

    print_abs_ratios <- function(y_in_raw, y_in, y_0_in_raw, y_0_in, i_L1_in, i_L2_in, i_L3_in, mask_in, skip) {
      print_template <- function(var_name, i_in, y_0_in_raw, y_0_in) {
        if(sum(i_in)==0) return(cat(paste0("\nNo examples of ", var_name, "!\n")))
        N <- length(i_in)
        cat(paste0("\nRatios for ", var_name, ":\n"))
        for(idx in c(1:3)) {
          a <- paste0(col_names[mask_in][[((idx*2)-1)]],"/",col_names[mask_in][[(idx*2)]])
          b <- sum(y_0_in[idx,i_in])/N
          c <- sum(y_0_in_raw[((idx*2)-1),i_in])/N
          d <- sum(y_0_in_raw[(idx*2),i_in])/N
          print_message(a, b, c, d)
        }
      }
      cat("\nDefault Ratios:\n")
      for(idx in c(1:3)) {
        a <- paste0(col_names[mask_in][[((idx*2)-1)]],"/",col_names[mask_in][[(idx*2)]])
        b <- y_in[idx]
        c <- y_in_raw[((idx*2)-1)]
        d <- y_in_raw[(idx*2)]
        print_message(a, b, c, d)
      }
      N <- length(i_L1_in)
      if(skip == 0) cat("\nRatios for ALL:\n")
      if(skip != 0) cat(paste0("\nRatios for ALL (All == DGA Status ", skip, "):\n"))
      for(idx in c(1:3)) {
        a <- paste0(col_names[mask_in][[((idx*2)-1)]],"/",col_names[mask_in][[(idx*2)]])
        b <- sum(y_0_in[idx,])/N
        c <- sum(y_0_in_raw[((idx*2)-1),])/N
        d <- sum(y_0_in_raw[(idx*2),])/N
        print_message(a, b, c, d)
      }
      if((verbose == TRUE) & (skip == 0)) {
        print_template("SL1", i_L1_in, y_0_in_raw, y_0_in)
        print_template("SL2", i_L2_in, y_0_in_raw, y_0_in)
        print_template("SL3", i_L3_in, y_0_in_raw, y_0_in)
      }
    }

    print_D <- function(d_name, D_L0, D_L1, D_L2, D_L3, D_Ln1, skip) {
      if(skip == 0) cat(paste0("\n(", d_name, "): Diagnoses for ALL:\n"))
      if(skip != 0) cat(paste0("\n(", d_name, "): Diagnoses for ALL (All == DGA Status ", skip, "):\n"))
      print(D_L0)
      if(skip != 0) return()
      if(is.nan(sum(D_L1))) {
        cat(paste0("\n(", d_name, "): No examples of SL1!\n"))
      } else {
        cat(paste0("\n(", d_name, "): Diagnoses for SL1:\n"))
        print(D_L1)
      }
      if(verbose == TRUE) {
        if(is.nan(sum(D_L2))) {
          cat(paste0("\n(", d_name, "): No examples of SL2!\n"))
        } else {
          cat(paste0("\n(", d_name, "): Diagnoses for SL2:\n"))
          print(D_L2)
        }
        if(is.nan(sum(D_L3))) {
          cat(paste0("\n(", d_name, "): No examples of SL3!\n"))
        } else {
          cat(paste0("\n(", d_name, "): Diagnoses for SL3:\n"))
          print(D_L3)
        }
      }
      if(is.nan(sum(D_Ln1))) {
        cat(paste0("\n(", d_name, "): No examples of SL2 or SL3!\n"))
      } else {
        cat(paste0("\n(", d_name, "): Diagnoses for SL2 and SL3:\n"))
        print(D_Ln1)
      }
    }
  }

  g <- dim(y_1)[1] # y_1 is [gas, samples]
  n <- dim(y_1)[2]
  gas_vals <- gas_data[1,-1]
  # Only plot sub-categories if there are multiple categories present
  skip <- ((sum(i_L1)==n))+((sum(i_L2)==n)*2)+((sum(i_L3)==n)*3)

  if((diagnosis == 1) | (diagnosis == 2)) { # DT
    # DT1 wants C2H4, CH4, C2H2
    mask <- match(c("C2H4", "CH4", "C2H2"), col_names)
    y   <- gas_vals[mask]
    y_0 <- prop.table(y_1[mask,], margin = 2)
    if(silent == FALSE) print_rel_ratios(y, y_0, mask, i_L1, i_L2, i_L3, skip)

    D_levels <- c("ND","PD","T1","T2","T3","DT","D1","D2")
    D1_0 <- factor(Interpret_DTn("1", y, 1), levels = D_levels)
    D1_0 <- table(D1_0)
    DT1_default <- prop.table(D1_0)
    if(silent == FALSE) {
      cat("\n(DT1): Default Diagnosis:\n")
      print(D1_0)
    }

    D1 <- factor(Interpret_DTn("1", y_0, n), levels = D_levels)
    DT1_L0  <- prop.table(table(D1))
    DT1_L1  <- prop.table(table(D1[i_L1]))
    DT1_L2  <- prop.table(table(D1[i_L2]))
    DT1_L3  <- prop.table(table(D1[i_L3]))
    DT1_Ln1 <- prop.table(table(D1[i_L2|i_L3]))
    if(silent == FALSE) print_D("DT1" , DT1_L0, DT1_L1, DT1_L2, DT1_L3, DT1_Ln1, skip)
    output$DT1_default_data <- y / sum(y)
    output$DT1_0_data <- y_0
    output$DT1_L0 <- DT1_L0
    output$DT1_L1 <- DT1_L1
    output$DT1_L2 <- DT1_L2
    output$DT1_L3 <- DT1_L3
    output$DT1_Ln1 <- DT1_Ln1

    # DT1 to DT4 and DT5 Referral Scenarios:
    # DT1 = {T3}      ->       DT5
    # DT1 = {T2}      -> DT4 & DT5
    # DT1 = {T1 | PD} -> DT4

    # DT4 wants CH4, H2, C2H6
    mask <- match(c("CH4", "H2", "C2H6"), col_names)
    y   <- gas_vals[mask]
    y_0 <- prop.table(y_1[mask,], margin = 2)
    if(silent == FALSE) print_rel_ratios(y, y_0, mask, i_L1, i_L2, i_L3, skip)

    D_levels <- c("NA","ND","PD","S","O","C")
    D <- factor(Interpret_DTn("4", y, 1), levels = D_levels)
    D[!(D1_0 %in% c("PD", "T1", "T2"))] <- factor("NA", levels = D_levels) # DT4 applies only when DT1 is T1, T2 or PD
    D <- table(D)
    DT4_default <- prop.table(D)
    if(silent == FALSE) {
      cat("\n(DT4): Default Diagnosis:\n")
      print(D)
    }

    D <- factor(Interpret_DTn("4", y_0, n), levels = D_levels)
    D[!(D1 %in% c("PD", "T1", "T2"))] <- factor("NA", levels = D_levels) # DT4 applies only when DT1 is T1, T2 or PD
    DT4_L0  <- prop.table(table(D))
    DT4_L1  <- prop.table(table(D[i_L1]))
    DT4_L2  <- prop.table(table(D[i_L2]))
    DT4_L3  <- prop.table(table(D[i_L3]))
    DT4_Ln1 <- prop.table(table(D[i_L2|i_L3]))
    if(silent == FALSE) print_D("DT4", DT4_L0, DT4_L1, DT4_L2, DT4_L3, DT4_Ln1, skip)
    output$DT4_default_data <- y / sum(y)
    output$DT4_0_data <- y_0
    output$DT4_L0 <- DT4_L0
    output$DT4_L1 <- DT4_L1
    output$DT4_L2 <- DT4_L2
    output$DT4_L3 <- DT4_L3
    output$DT4_Ln1 <- DT4_Ln1

    # DT5 wants C2H4, CH4, C2H6
    mask <- match(c("C2H4", "CH4", "C2H6"), col_names)
    y   <- gas_vals[mask]
    y_0 <- prop.table(y_1[mask,], margin = 2)
    if(silent == FALSE) print_rel_ratios(y, y_0, mask, i_L1, i_L2, i_L3, skip)

    D_levels <- c("NA","ND","PD","S","O","C","T2","T3")
    D <- factor(Interpret_DTn("5", y, 1), levels = D_levels)
    D[!(D1_0 %in% c("T2", "T3"))] <- factor("NA", levels = D_levels) # DT5 applies only when DT1 is T2 or T3
    D <- table(D)
    DT5_default <- prop.table(D)
    if(silent == FALSE) {
      cat("\n(DT5): Default Diagnosis:\n")
      print(D)
    }

    D <- factor(Interpret_DTn("5", y_0, n), levels = D_levels)
    D[!(D1 %in% c("T2", "T3"))] <- factor("NA", levels = D_levels) # DT5 applies only when DT1 is T2 or T3
    DT5_L0  <- prop.table(table(D))
    DT5_L1  <- prop.table(table(D[i_L1]))
    DT5_L2  <- prop.table(table(D[i_L2]))
    DT5_L3  <- prop.table(table(D[i_L3]))
    DT5_Ln1 <- prop.table(table(D[i_L2|i_L3]))
    if(silent == FALSE) print_D("DT5", DT5_L0, DT5_L1, DT5_L2, DT5_L3, DT5_Ln1, skip)
    output$DT5_default_data <- y / sum(y)
    output$DT5_0_data <- y_0
    output$DT5_L0 <- DT5_L0
    output$DT5_L1 <- DT5_L1
    output$DT5_L2 <- DT5_L2
    output$DT5_L3 <- DT5_L3
    output$DT5_Ln1 <- DT5_Ln1
  }

  if((diagnosis == 1) | (diagnosis == 3)) { # RR
    # RR wants C2H2:C2H4, CH4:H2, C2H4:C2H6
    mask <- match(c("C2H2", "C2H4", "CH4", "H2", "C2H4", "C2H6"), col_names)
    y <- unlist(gas_vals[mask])
    y[y == 0] <- 0.00001 # To avoid /0 error
    y <- c(y[1]/y[2], y[3]/y[4], y[5]/y[6])
    names(y) <- c("C2H2_C2H4", "CH4_H2", "C2H4_C2H6")

    y_0 <- y_1[mask,]
    y_0[y_0 == 0] <- 0.00001 # To avoid /0 error
    y_0 <- t(array(c(y_0[1,]/y_0[2,], y_0[3,]/y_0[4,], y_0[5,]/y_0[6,]), dim = c(n,3)))

    if((silent == FALSE)) print_abs_ratios(unlist(gas_vals[mask]), y, y_1[mask,], y_0, i_L1, i_L2, i_L3, mask, skip)

    D_levels <- c("ND","C0","C1","C2","C3","C4","C5")
    D <- factor(Interpret_RR(y, 1), levels = D_levels)
    if((silent == FALSE)) output$RR_D_default <- D # Just for RR plotting, want actual diagnosis for each
    D <- table(D)
    RR_default <- prop.table(D)
    if(silent == FALSE) {
      cat("\nRR: Default Diagnosis:\n")
      print(D)
    }

    D <- factor(Interpret_RR(y_0, n), levels = D_levels)
    if((silent == FALSE)) output$RR_D <- D # Just for RR plotting, want actual diagnosis for each
    RR_L0  <- prop.table(table(D))
    RR_L1  <- prop.table(table(D[i_L1]))
    RR_L2  <- prop.table(table(D[i_L2]))
    RR_L3  <- prop.table(table(D[i_L3]))
    RR_Ln1 <- prop.table(table(D[i_L2|i_L3]))
    if(silent == FALSE) print_D("RR" , RR_L0, RR_L1, RR_L2, RR_L3, RR_Ln1, skip)
    output$Abs_default_data <- y
    output$Abs_0_data <- y_0
    output$RR_L0 <- RR_L0
    output$RR_L1 <- RR_L1
    output$RR_L2 <- RR_L2
    output$RR_L3 <- RR_L3
    output$RR_Ln1 <- RR_Ln1
  }

  if((diagnosis == 1) | (diagnosis == 4)) { # IEC
    # IEC wants C2H2:C2H4, CH4:H2, C2H4:C2H6
    if(diagnosis == 4) { # Will already be available from RR otherwise
      mask <- match(c("C2H2", "C2H4", "CH4", "H2", "C2H4", "C2H6"), col_names)
      y <- unlist(gas_vals[mask])
      y[y == 0] <- 0.00001 # To avoid /0 error
      y <- c(y[1]/y[2], y[3]/y[4], y[5]/y[6])
      names(y) <- c("C2H2_C2H4", "CH4_H2", "C2H4_C2H6")

      y_0 <- y_1[mask,]
      y_0[y_0 == 0] <- 0.00001 # To avoid /0 error
      y_0 <- t(array(c(y_0[1,]/y_0[2,], y_0[3,]/y_0[4,], y_0[5,]/y_0[6,]), dim = c(n,3)))

      output$Abs_default_data <- y
      output$Abs_0_data <- y_0

      if((silent == FALSE)) print_abs_ratios(unlist(gas_vals[mask]), y, y_1[mask,], y_0, i_L1, i_L2, i_L3, mask, skip)
    }

    # IEC 1
    D_levels <- c("ND","PD","D1","D2","T1","T2","T3")
    D <- factor(Interpret_IECn("1", y, 1), levels = D_levels)
    if((silent == FALSE)) output$IEC1_D_default <- D # Just for IEC plotting, want actual diagnosis for each
    D <- table(D)
    IEC1_default <- prop.table(D)
    if(silent == FALSE) {
      cat("\nIEC 1: Default Diagnosis:\n")
      print(D)
    }

    D <- factor(Interpret_IECn("1", y_0, n), levels = D_levels)
    if((silent == FALSE)) output$IEC1_D <- D # Just for IEC plotting, want actual diagnosis for each
    IEC1_L0  <- prop.table(table(D))
    IEC1_L1  <- prop.table(table(D[i_L1]))
    IEC1_L2  <- prop.table(table(D[i_L2]))
    IEC1_L3  <- prop.table(table(D[i_L3]))
    IEC1_Ln1 <- prop.table(table(D[i_L2|i_L3]))
    if(silent == FALSE) print_D("IEC 1", IEC1_L0, IEC1_L1, IEC1_L2, IEC1_L3, IEC1_Ln1, skip)
    output$IEC1_L0 <- IEC1_L0
    output$IEC1_L1 <- IEC1_L1
    output$IEC1_L2 <- IEC1_L2
    output$IEC1_L3 <- IEC1_L3
    output$IEC1_Ln1 <- IEC1_Ln1

    # IEC 2
    D_levels <- c("ND","PD","D","T")
    D <- factor(Interpret_IECn("2", y, 1), levels = D_levels)
    if((silent == FALSE)) output$IEC2_D_default <- D # Just for IEC plotting, want actual diagnosis for each
    D <- table(D)
    IEC2_default <- prop.table(D)
    if(silent == FALSE) {
      cat("\nIEC 2: Default Diagnosis:\n")
      print(D)
    }

    D <- factor(Interpret_IECn("2", y_0, n), levels = D_levels)
    if((silent == FALSE)) output$IEC2_D <- D # Just for IEC plotting, want actual diagnosis for each
    IEC2_L0  <- prop.table(table(D))
    IEC2_L1  <- prop.table(table(D[i_L1]))
    IEC2_L2  <- prop.table(table(D[i_L2]))
    IEC2_L3  <- prop.table(table(D[i_L3]))
    IEC2_Ln1 <- prop.table(table(D[i_L2|i_L3]))
    if(silent == FALSE) print_D("IEC 2", IEC2_L0, IEC2_L1, IEC2_L2, IEC2_L3, IEC2_Ln1, skip)
    output$IEC2_L0 <- IEC2_L0
    output$IEC2_L1 <- IEC2_L1
    output$IEC2_L2 <- IEC2_L2
    output$IEC2_L3 <- IEC2_L3
    output$IEC2_Ln1 <- IEC2_Ln1
  }

  return(output)
}
