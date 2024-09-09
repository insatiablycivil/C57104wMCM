#' Define_RR: Defines zones for plotting Rogers Ratio Diagnostic Logic as per IEEE Table 5
#'
#' This defines zones for plotting Rogers Ratio in accordance to IEEE C57.104-2019 Table 5.
#'
#' @return Data.frame:
#' \itemize{
#'  \item \code{xmin} | \code{xmax} | \code{ymin} | \code{ymax}: Numeric. Used to define rectangles in a ggplot.
#'  \item \code{label}: String. Used for colour scaling in ggplot.
#'  \item \code{name}: String. Used for defining which facet to plot in.
#'  }
#' @export
Define_RR <- function() {
  # These define the diagnostic zones for plotting
  # These create rectangular shaded zones on a facet_grid plot of 2 plots.
  zones <- data.frame(
    #         Truncated         Unknown      C0   C1   C2   C3   C4   C5
    xmin = c( 0 ,10 , 0  , 0  , 0  ,0  ,NA , 0  , 0  , 3  , 1  , 1  , 3  ,
              0 ,10 , 0  , 0  , 0  ,0  ,1  , 0.1, 0  , 0.1, 0.1, 1  , 1  ),
    xmax = c(10 ,Inf,10  ,.005, Inf,3  ,NA , 1  , 1  , Inf, 3  , 3  , Inf,
             10 ,Inf,10  ,.005, Inf,0.1,Inf, 1  , 0.1, 1  , 1  , Inf, Inf),
    ymin = c(10 , 0 , 0  ,.005, 3  ,0.1,NA , 0  , 0  , 0.1, 0  , 0  , 0  ,
             10 , 0 , 0  ,.005, 3  ,0.1,0.1, 0  , 0  , 0.1, 0  , 0  , 0  ),
    ymax = c(Inf,Inf,.005,10  , Inf,3  ,NA , 0.1, 0.1, 3.0, 0.1, 0.1, 0.1,
             Inf,Inf,.005,10  , Inf,3  ,3  , 0.1, 0.1, 3.0, 0.1, 0.1, 0.1),
    label = rep(c(rep('Truncated', 4), rep('.Unknown', 3), 'C0: Normal',
                  'C1: PD','C2: Arcing','C3: Thermal Low',
                  'C4: Thermal<700C','C5: Thermal>700C'),2),
    name = c(rep("C2H4/C2H6",13), rep("CH4/H2",13)))
  return(zones)
}
