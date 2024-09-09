#' Define_IECn: Defines zones for plotting Duval Triangles 1, 4, or 5 Diagnostic Logic as per IEEE Table 6, D.3, and D.4, respectively.
#'
#' This defines zones for plotting in accordance to either IEEE C57.104-2019 Table 6 or Table D.3, or Table D.4. Table D.3 is modified here such that original definition for C: (CH4 >= 36 & C2H6 >= 24) is instead: (CH4 >= 36 & C2H6 < 24).
#'
#' @param DT_n String: Value corresponding to either Triangle 1, 4, or 5. Expected value: 1, 4, or 5.
#' @return List:
#' \itemize{
#'  \item \code{labels}: String. Used to define axes used in ternary plot.
#'  \item \code{zones}: Numeric. Used for defining diagnostic zones used in ternary plot.
#'  \item \code{colours}: RGB. Used for defining colours of zones.
#'  }
#' @export
Define_DTn <- function(DT_n) {
  # These define the diagnostic zones for plotting
  Define_DT1 <- function() {
    # Encoded fault zones
    # PD - Corona  LHS  RHS  BOT
    PD <- matrix(c(100, 000, 000,
                   098, 000, 002,
                   098, 002, 000),
                 ncol = 3, byrow = TRUE)
    # T1 - Thermal Faults of Temp < 300C
    T1 <- matrix(c(098, 000, 002,
                   096, 000, 004,
                   076, 020, 004,
                   080, 020, 000,
                   098, 002, 000),
                 ncol = 3, byrow = TRUE)
    # T2 - Thermal Faults of 300C < Temp < 700C
    T2 <- matrix(c(080, 020, 000,
                   076, 020, 004,
                   046, 050, 004,
                   050, 050, 000),
                 ncol = 3, byrow = TRUE)
    # T3 - Thermal Faults of Temp > 700C
    T3 <- matrix(c(050, 050, 000,
                   035, 050, 015,
                   000, 085, 015,
                   000, 100, 000),
                 ncol = 3, byrow = TRUE)
    # DT - Mixtures of Electrical and Thermal Faults
    DT <- matrix(c(096, 000, 004,
                   087, 000, 013,
                   047, 040, 013,
                   031, 040, 029,
                   000, 071, 029,
                   000, 085, 015,
                   035, 050, 015,
                   046, 050, 004),
                 ncol = 3, byrow = TRUE)
    # D1 - Electrical Dischanrge of Low Energy
    D1 <- matrix(c(087, 000, 013,
                   000, 000, 100,
                   000, 023, 077,
                   064, 023, 013),
                 ncol = 3, byrow = TRUE)
    # D2 - Electrical Discharge of High Energy
    D2 <- matrix(c(064, 023, 013,
                   000, 023, 077,
                   000, 071, 029,
                   031, 040, 029,
                   047, 040, 013),
                 ncol = 3, byrow = TRUE)

    labels <- list(point = 'up',
                   atip = 'CH4',
                   btip = 'C2H4',
                   ctip = 'C2H2',
                   alab = '% CH4',
                   blab = '% C2H4',
                   clab = '% C2H2')
    zones <- list(PD = PD,
                  T1 = T1,
                  T2 = T2,
                  T3 = T3,
                  DT = DT,
                  D1 = D1,
                  D2 = D2)
    colours <- list(PD = rgb(1,1,1,0.5),
                    T1 = rgb(0.89803,0.92157,0.11765,0.5),
                    T2 = rgb(0.89412,0.66667,0.13725,0.5),
                    T3 = rgb(0.96078,0.43922,0.12549,0.5),
                    DT = rgb(0.75294,0.58039,0.33725,0.5),
                    D1 = rgb(0.74902,0.26275,0.14902,0.5),
                    D2 = rgb(1,0,0,0.5))

    return(list(labels = labels, zones = zones, colours = colours))
  }

  Define_DT4 <- function() {
    # Encoded fault zones
    # PD - Corona       LHS  RHS  BOT
    PD <- matrix(c(098, 002, 000,
                   097, 002, 001,
                   084, 015, 001,
                   085, 015, 000),
                 ncol = 3, byrow = TRUE)
    # S - Stray Gassing of Mineral Oil (Temp < 200C)
    S <-  matrix(c(100, 000, 000,
                   054, 000, 046,
                   009, 045, 046,
                   009, 061, 030,
                   015, 055, 030,
                   015, 061, 024,
                   040, 036, 024,
                   064, 036, 000,
                   085, 015, 000,
                   084, 015, 001,
                   097, 002, 001,
                   098, 002, 000),
                 ncol = 3, byrow = TRUE)
    # O - Overheating (Temp < 250C)
    O <-  matrix(c(009, 000, 091,
                   000, 000, 100,
                   000, 070, 030,
                   009, 061, 030),
                 ncol = 3, byrow = TRUE)
    # C - Hot Spots w Carbonisation of Paper (Temp > 300C)
    C <-  matrix(c(064, 036, 000,
                   040, 036, 024,
                   015, 061, 024,
                   015, 055, 030,
                   000, 070, 030,
                   000, 100, 000),
                 ncol = 3, byrow = TRUE)
    # ND - Not Determined
    ND <- matrix(c(054, 000, 046,
                   009, 000, 091,
                   009, 045, 046),
                 ncol = 3, byrow = TRUE)

    labels <- list(point = 'up',
                   atip = 'H2',
                   btip = 'CH4',
                   ctip = 'C2H6',
                   alab = '% H2',
                   blab = '% CH4',
                   clab = '% C2H6')
    zones <- list(PD = PD,
                  S = S,
                  O = O,
                  C = C,
                  ND = ND)
    colours <- list(PD = rgb(1,1,1,0.5),
                    S =  rgb(0.89803,0.92157,0.11765,0.5),
                    O =  rgb(0.75000,0.82000,0.07000,0.5),
                    C =  rgb(0.75294,0.58039,0.33725,0.5),
                    ND = rgb(0.74910,0.74910,0.74910,0.5))

    return(list(labels = labels, zones = zones, colours = colours))
  }

  Define_DT5 <- function() {
    # Encoded fault zones
    # PD - Corona  LHS  RHS  BOT
    PD <- matrix(c(098, 000, 002,
                   086, 000, 014,
                   085, 001, 014,
                   097, 001, 002),
                 ncol = 3, byrow = TRUE)
    # O - Overheating (Temp < 250C) # TWO REGIONS
    O_1 <-  matrix(c(100, 000, 000,
                     098, 000, 002,
                     097, 001, 002,
                     085, 001, 014,
                     076, 010, 014,
                     090, 010, 000),
                   ncol = 3, byrow = TRUE)
    # O - Overheating (Temp < 250C) # TWO REGIONS
    O_2 <-  matrix(c(046, 000, 054,
                     000, 000, 100,
                     000, 010, 090,
                     036, 010, 054),
                   ncol = 3, byrow = TRUE)
    # S - Stray Gassing of Mineral Oil (Temp < 200C)
    S <-  matrix(c(086, 000, 014,
                   046, 000, 054,
                   036, 010, 054,
                   076, 010, 014),
                 ncol = 3, byrow = TRUE)
    # T2 - Thermal Faults of 300C < Temp < 700C
    T2 <-  matrix(c(090, 010, 000,
                    078, 010, 012,
                    053, 035, 012,
                    065, 035, 000),
                  ncol = 3, byrow = TRUE)
    # T3 - Thermal Faults of Temp > 700C
    T3 <-  matrix(c(065, 035, 000,
                    053, 035, 012,
                    038, 050, 012,
                    036, 050, 014,
                    016, 070, 014,
                    000, 070, 030,
                    035, 035, 030,
                    000, 035, 065,
                    000, 100, 000),
                  ncol = 3, byrow = TRUE)
    # C - Hot Spots w Carbonisation of Paper (Temp > 300C)
    C <-  matrix(c(078, 010, 012,
                   060, 010, 030,
                   000, 070, 030,
                   016, 070, 014,
                   036, 050, 014,
                   038, 050, 012),
                 ncol = 3, byrow = TRUE)
    # ND - Not Determined
    ND <- matrix(c(060, 010, 030,
                   000, 010, 090,
                   000, 035, 065,
                   035, 035, 030),
                 ncol = 3, byrow = TRUE)

    labels <- list(point='up',
                   atip='CH4',
                   btip='C2H4',
                   ctip='C2H6',
                   alab='% CH4',
                   blab='% C2H4',
                   clab='% C2H6')
    zones <- list(PD = PD,
                  O_1 = O_1,
                  O_2 = O_2,
                  S = S,
                  T2 = T2,
                  T3 = T3,
                  C = C,
                  ND = ND)
    colours <- list(PD =  rgb(1,1,1,0.5),
                    O_1 = rgb(0.75000,0.82000,0.07000,0.5),
                    O_2 = rgb(0.75000,0.82000,0.07000,0.5),
                    S =   rgb(0.89803,0.92157,0.11765,0.5),
                    T2 =  rgb(0.89412,0.66667,0.13725,0.5),
                    T3 =  rgb(0.96078,0.43922,0.12549,0.5),
                    C =   rgb(0.75294,0.58039,0.33725,0.5),
                    ND =  rgb(0.74910,0.74910,0.74910,0.5))

    return(list(labels = labels, zones = zones, colours = colours))
  }

  if(DT_n == "1") return(Define_DT1())
  if(DT_n == "4") return(Define_DT4())
  if(DT_n == "5") return(Define_DT5())
}
