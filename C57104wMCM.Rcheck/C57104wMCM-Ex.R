pkgname <- "C57104wMCM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('C57104wMCM')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("C57104wMCM")
### * C57104wMCM

flush(stderr()); flush(stdout())

### Name: C57104wMCM
### Title: C57104wMCM: Simplified C57.104 DGA Screening with propagated
###   uncertainty via a Monte Carlo Method
### Aliases: C57104wMCM

### ** Examples

# Example 1: Two gases, two samples, individually assigned accuracies, no correlation, and no diagnostics
# gas_data: Dates and gas values (PPM). With only two samples, Table 4 will not be used.
gas_data_in <- data.frame(date = as.Date(c("2024-06-20","2024-06-24")),
                          H2 = c(8,4), CH4 = c(80,160))
# limit_data: IEEE Tables 1-4 Limits for each gas (PPM and PPM/Year).
limit_data_in <- data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03))
# accuracy: Multiplied by gas value to give uncertainty.
# If "Normal", this represents 95% probability level, else 100%.
# Can also give a single value to be used for all gases
accuracy_in <- setNames(c(0.15, 0.50), c("H2", "CH4"))
# min_accuracy: Treated as PPM and used if the above would give an otherwise smaller value.
# Either can be given as NA to force the use of the other. This too can be given as single value.
min_accuracy_in = setNames(c(1, 5), c("H2", "CH4"))
C57104wMCM(gas_data = gas_data_in,
           limit_data = limit_data_in,
           accuracy = accuracy_in,
           min_accuracy = min_accuracy_in,
           distribution_shape = "N",
           rho = 0.5,
           diagnosis = 0,
           n = 1000,
           flex = FALSE, silent = FALSE, verbose = FALSE, validate = FALSE)

# Example 2: Seven gases, five samples, shared accuracies, specified correlation matrix, and Duval Triangle diagnostics
gas_data_in <- data.frame(date=as.Date(c("2024-06-20","2024-06-21","2024-06-22","2024-06-23","2024-06-24")),
                          H2=c(80,80,80,80,82), CH4=c(90,95,100,105,107), C2H6=c(80,100,80,80,82),
                          C2H4=c(30,30,31,29,31), C2H2=c(1,1,0,1,3), CO=c(800,850,860,850,852),
                          CO2=c(900,950,940,940,942))
limit_data_in <- data.frame(H2=c(80,200,40,20),CH4=c(90,150,30,10),C2H6=c(90,175,25,9),
                            C2H4=c(50,100,20,7),C2H2=c(1,2,0.5,0.5),CO=c(900,1100,20,100),
                            CO2=c(9000,12500,250,1000))
# rho: If given as a matrix, it is interpreted as a correlation matrix, else as rho.
# In this example, the correlation matrix is not positive definite.
# However, the flex argument allows for it to be adjusted to make it so.
rho_in <- matrix(rep(1,49), ncol=7)
C57104wMCM(gas_data = gas_data_in,
           limit_data = limit_data_in,
           accuracy = 0.15,
           min_accuracy = 1,
           distribution_shape = "T",
           rho = rho_in,
           diagnosis = "DT",
           n = 500,
           flex = TRUE, silent = FALSE, verbose = TRUE, validate = TRUE)



cleanEx()
nameEx("Diagnosis")
### * Diagnosis

flush(stderr()); flush(stdout())

### Name: Diagnosis
### Title: Diagnosis: Apply Diagnostics in accordance to either Duval
###   Triangles, Rogers Ratios, or IEC Tables, or all.
### Aliases: Diagnosis

### ** Examples

# y_1: Numeric array of generated samples, shaped: Gases, Samples
y_1_in <- array(c(82.39, 107.51,  82.39,  31.15,   3.03, 856.02, 946.45,  77.95, 101.71,  77.95,  29.47,   2.67, 809.89, 895.44,
                  78.98, 103.06,  78.98,  29.86,   2.75, 820.60, 907.28,  84.62, 110.42,  84.62,  31.99,   3.21, 879.26, 972.14,
                  74.19,  96.81,  74.19,  28.05,   2.37, 770.88, 852.31,  85.56, 111.64,  85.56,  32.35,   3.29, 888.97, 982.88,
                  84.71, 110.54,  84.71,  32.03,   3.22, 880.20, 973.18,  74.09,  96.68,  74.09,  28.01,   2.36, 769.80, 851.12,
                  83.24, 108.62,  83.24,  31.47,   3.10, 864.91, 956.27,  76.22,  99.45,  76.22,  28.81,   2.53, 791.90, 875.55), dim = c(7,10))
gas_data_in <- data.frame(date=c(4,3,2,1,0),
                          H2=c(80,80,80,80,82), CH4=c(90,95,100,105,107), C2H6=c(80,100,80,80,82),
                          C2H4=c(30,30,31,29,31), C2H2=c(1,1,0,1,3), CO=c(800,850,860,850,852),
                          CO2=c(900,950,940,940,942))
# temp_rndm_screen is just for generating screening mask as an example.
temp_rndm_screen <- sample(1:3, 10, replace=TRUE)
i_L1_in <- temp_rndm_screen == 1
i_L2_in <- temp_rndm_screen == 2
i_L3_in <- temp_rndm_screen == 3

Diagnosis(y_1 = y_1_in,
          gas_data = gas_data_in,
          diagnosis = 1,
          col_names = c("H2", "CH4", "C2H6", "C2H4", "C2H2", "CO", "CO2"),
          i_L1 = i_L1_in,
          i_L2 = i_L2_in,
          i_L3 = i_L3_in,
          silent = FALSE, verbose = FALSE)



cleanEx()
nameEx("Diagnosis_Plot")
### * Diagnosis_Plot

flush(stderr()); flush(stdout())

### Name: Diagnosis_Plot
### Title: Diagnosis_Plot: Plotting outputs from Diagnosis function for
###   validation.
### Aliases: Diagnosis_Plot

### ** Examples

# Not recreating required inputs to run this example
# Diagnosis_Plot(DI = diagnosed,
#                SC = screened,
#                PP = preprocessed)



cleanEx()
nameEx("Interpret_DTn")
### * Interpret_DTn

flush(stderr()); flush(stdout())

### Name: Interpret_DTn
### Title: Interpret_DTn: Applies Duval Triangles 1, 4, or 5 Diagnostic
###   Logic as per IEEE Table 6, D.3, and D.4, respectively.
### Aliases: Interpret_DTn

### ** Examples

Interpret_DTn(DT_n = "4",
               y_in = array(c(0.090, 0.002, 0.401, 0.120, 0.001, 0.310), dim = c(3,2)),
               n_in = 2)



cleanEx()
nameEx("Interpret_IECn")
### * Interpret_IECn

flush(stderr()); flush(stdout())

### Name: Interpret_IECn
### Title: Interpret_IECn: Applies IEC Table 1 or Table 2 Diagnostic Logic
### Aliases: Interpret_IECn

### ** Examples

Interpret_IECn(IEC_n = "1",
               y_in = array(c(0.090, 0.002, 0.401, 0.120, 0.001, 0.310), dim = c(3,2)),
               n_in = 2)



cleanEx()
nameEx("Interpret_RR")
### * Interpret_RR

flush(stderr()); flush(stdout())

### Name: Interpret_RR
### Title: Interpret_RR: Applies Rogers Ratio Diagnostic Logic as per IEEE
###   Table 5
### Aliases: Interpret_RR

### ** Examples

Interpret_RR(y_in = array(c(0.090, 0.002, 0.401, 0.120, 0.001, 0.310), dim = c(3,2)),
             n_in = 2)



cleanEx()
nameEx("Preprocess")
### * Preprocess

flush(stderr()); flush(stdout())

### Name: Preprocess
### Title: Preprocess: Validate inputs and preprocess them into expected
###   data formats / values
### Aliases: Preprocess

### ** Examples

Preprocess(gas_data = data.frame(date = as.Date(c("2024-06-20","2024-06-21","2024-06-22","2024-06-23","2024-06-24")), H2 = c(8,4,2,1,0.5), CH4 = c(10,20,40,80,160)),
           limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
           accuracy = setNames(c(0.15, 0.50), c("H2", "CH4")),
           min_accuracy = setNames(c(1, 5), c("H2", "CH4")),
           distribution_shape = "N",
           rho = 0.5,
           diagnosis = 0,
           n = 1000,
           flex = FALSE, silent = FALSE, verbose = FALSE, validate = FALSE)
Preprocess(gas_data = data.frame(date=as.Date(c("2024-06-20","2024-06-21","2024-06-22","2024-06-23","2024-06-24")), H2=c(80,80,80,80,82),CH4=c(90,95,100,105,107),C2H6=c(80,100,80,80,82),C2H4=c(30,30,31,29,31),C2H2=c(1,1,0,1,3),CO=c(800,850,860,850,852),CO2=c(900,950,940,940,942)),
           limit_data = data.frame(H2=c(80,200,40,20),CH4=c(90,150,30,10),C2H6=c(90,175,25,9),C2H4=c(50,100,20,7),C2H2=c(1,2,0.5,0.5),CO=c(900,1100,20,100),CO2=c(9000,12500,250,1000)),
           accuracy = 0.15,
           min_accuracy = 1,
           distribution_shape = "T",
           rho = matrix(rep(1,49), ncol=7),
           diagnosis = "DT",
           n = 1000,
           flex = TRUE, silent = FALSE, verbose = TRUE, validate = TRUE)



cleanEx()
nameEx("Sample")
### * Sample

flush(stderr()); flush(stdout())

### Name: Sample
### Title: Sample: Randomly generate samples to later apply Screening to
###   for Monte Carlo estimation
### Aliases: Sample

### ** Examples

Sample(gas_data = data.frame(date = c(4,3,2,1,0), H2 = c(8,4,2,1,0.5), CH4 = c(10,20,40,80,160)),
       limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
       accuracy = setNames(c(0.15, 0.50), c("H2", "CH4")),
       min_accuracy = setNames(c(1, 5), c("H2", "CH4")),
       distribution_shape = 1,
       rho = 0.5,
       cor_mat = NA,
       n = 1000,
       col_names = c("H2", "CH4"),
       sd_data = data.frame(date = c(4,3,2,1,0), H2 = c(1.2,0.6,0.3,0.15,0.075), CH4 = c(5,10,20,40,80)))
Sample(gas_data = data.frame(date = c(4,3,2,1,0), H2 = c(8,4,2,1,0.5), CH4 = c(10,20,40,80,160)),
       limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
       accuracy = 0.15,
       min_accuracy = 1,
       distribution_shape = 1,
       rho = NA,
       cor_mat = matrix(c(1,0,0,1), ncol=2),
       n = 1000,
       col_names = c("H2", "CH4"),
       sd_data = data.frame(date = c(4,3,2,1,0), H2 = c(1.2,0.6,0.3,0.15,0.075), CH4 = c(0.75,1.50,3,6,12)))



cleanEx()
nameEx("Sample_Plot")
### * Sample_Plot

flush(stderr()); flush(stdout())

### Name: Sample_Plot
### Title: Sample_Plot: Plotting outputs from Sample function for
###   validation.
### Aliases: Sample_Plot

### ** Examples

# Not recreating required inputs to run this example
# Sample_Plot(y_1 = y_N[1,,],
#             PP = preprocessed)



cleanEx()
nameEx("Screening")
### * Screening

flush(stderr()); flush(stdout())

### Name: Screening
### Title: Screening: Apply Simplified IEEE C57.104-2019 Screening - Use
###   all samples for Tables 1-4.
### Aliases: Screening

### ** Examples

Screening(y_N = array(c(2,1,0.5,40,80,160,2.1,1.1,0.6,40.1,80.1,160.1), dim = c(2,2,2)),
          gas_data = data.frame(date = as.Date(c(2,1,0)), H2 = c(8,4,1), CH4 = c(40,80,160)),
          limit_data = data.frame(H2 = c(2, 3, 2, 0.05), CH4 = c(60,180,70,0.03)),
          silent = FALSE, verbose = FALSE)



cleanEx()
nameEx("Screening_Plot")
### * Screening_Plot

flush(stderr()); flush(stdout())

### Name: Screening_Plot
### Title: Screening_Plot: Plotting outputs from Screening function for
###   validation.
### Aliases: Screening_Plot

### ** Examples

# Not recreating required inputs to run this example
# Screening_Plot(SC = screened,
#                PP = preprocessed)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
