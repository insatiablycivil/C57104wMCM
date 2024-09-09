#' Define_IECn: Defines zones for plotting IEC Table 1 or Table 2 Diagnostic Logic
#'
#' This defines zones for plotting in accordance to either IEC 60599:2022 Table 1 or Table 2. For Table 1, bounds shown in Figure B.1 are used.
#'
#' @param IEC_n String: Value corresponding to either Table 1 or Table 2. Expected value: 1 or 2.
#' @return Data.frame:
#' \itemize{
#'  \item \code{xmin} | \code{xmax} | \code{ymin} | \code{ymax}: Numeric. Used to define rectangles in a ggplot.
#'  \item \code{label}: String. Used for colour scaling in ggplot.
#'  \item \code{name}: String. Used for defining which facet to plot in.
#'  }
#' @export
Define_IECn <- function(IEC_n) {
  # These define the diagnostic zones for plotting
  Define_IEC1 <- function() {
    # These create rectangular shaded zones on a facet_grid plot of 2 plots.
    zones <- data.frame(
      #        Truncated        Unknown          PD   Dn       Tn
      xmin = c(0  ,10 ,0  ,0  , 0  ,1  ,2  ,4  , 0  , 1  ,2  , 0  ,1  ,4  ,
               0  ,10 ,0  ,0  , 0  ,0.1,0.5,1  , 0  , 0.1,0.1, 1  ,1  ,1  ),
      xmax = c(10 ,Inf,10 ,.005,1  ,2  ,4  ,Inf, 0.2, Inf,Inf, 1  ,4  ,Inf,
               10 ,Inf,10 ,.005,0.1,1  ,1  ,Inf, 0.1, 0.5,  1, Inf,Inf,Inf),
      ymin = c(10 ,0  ,0  ,.005,.01,0.1,0.1,0.2, 0  , 1  ,0.6, 0  ,0  ,0  ,
               10 ,0  ,0  ,.005,.01,0  ,2.5,0.2, 0  , 1  ,0.6, 0  ,0  ,0  ),
      ymax = c(Inf,Inf,.005,10 ,Inf,1  ,0.6,0.6,0.01,Inf ,2.5,0.01,0.1,0.2,
               Inf,Inf,.005,10 ,Inf,0.6,Inf,Inf,0.01,Inf ,2.5,0.01,0.1,0.2),
      label = rep(c(rep('Truncated', 4), rep('.Unknown', 4), '1.PD',
                    '2.Arcing Low','3.Arcing High',
                    '4.Thermal<300C','5.Thermal<700C','6.Thermal>700C'),2),
      name = c(rep("C2H4/C2H6",14),rep("CH4/H2",14)))
    return(zones)
  }

  Define_IEC2 <- function() {
    # These create rectangular shaded zones. Only 1 plot needed.
    zones <- data.frame(
      #        Truncated         PD   D    T    NA
      xmin = c(0  ,10 ,0  ,0   , 0  , 0  , 0  , NA),
      xmax = c(10 ,Inf,10 ,.005, 0.2, Inf, Inf, NA),
      ymin = c(10 ,0  ,0  ,.005, 0  , 0.2, 0  , NA),
      ymax = c(Inf,Inf,.005,10 , Inf, Inf, 0.2, NA),
      label = c(rep('Truncated', 4), '1.PD','2.Arcing','3.Thermal','.NA'),
      name = rep("CH4/H2",8))
    return(zones)
  }

  if(IEC_n == "1") return(Define_IEC1())
  if(IEC_n == "2") return(Define_IEC2())
}
