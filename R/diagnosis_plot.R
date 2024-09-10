#' Diagnosis_Plot: Plotting outputs from Diagnosis function for validation.
#'
#' This is internal function to plot generated diagnosis outputs for validation.
#'
#' @param DI List: Outputs from Diagnosis function.
#' @param SC List: Outputs from Screening function.
#' @param PP List: Outputs from Preprocess function.
#' @examples
#' # Not recreating required inputs to run this example
#' # Diagnosis_Plot(DI = diagnosed,
#' #                SC = screened,
#' #                PP = preprocessed)
#' @export
Diagnosis_Plot <- function(DI, SC, PP) {
  # For viewing the data for validation purposes.
  # diagnosis_plot_<DTn|IECn|RR> contain functions demarcating diagnostic zones
  # Plots points graphically, and also as a frequency plot

  d <- dim(SC$T12_vals)[1]
  g <- dim(SC$T12_vals)[2]
  n <- dim(SC$T12_vals)[3]

  # Only plot sub-categories if there are multiple categories present
  skip <- ((sum(SC$L1_mask)==n))+((sum(SC$L2_mask)==n)*2)+((sum(SC$L3_mask)==n)*3)

  if((PP$diagnosis == 1) | (PP$diagnosis == 2)) { # DT 1|4|5
    library(Ternary)

    Summarise_DTn <- function(DT_n, samples, default_val, L1_mask, L2_mask, L3_mask, skip) {
      DTn <- Define_DTn(DT_n)
      Plot_DTn(DTn, DT_n, t(samples)*100, default_val*100, "0", skip)
      if(skip == 0) {
        if(sum(L1_mask)>0) Plot_DTn(DTn, DT_n, t(samples[,L1_mask])*100, default_val*100, "1", 0)
        # if(sum(L2_mask)>0) Plot_DTn(DTn, DT_n, t(samples[L2_mask])*100, default_val*100, "2", 0)
        # if(sum(L3_mask)>0) Plot_DTn(DTn, DT_n, t(samples[L3_mask])*100, default_val*100, "3", 0)
        if(sum(!L1_mask)>0) Plot_DTn(DTn, DT_n, t(samples[,!L1_mask])*100, default_val*100, "n1", 0)
      }
    }

    Plot_DTn <- function(DTn, DT_n, sampled, default, L_n, skip) {
      if(L_n == "0") {
        title_name <- switch(skip,
                             "All DGA Status (All == 1)",
                             "All DGA Status (All == 2)",
                             "All DGA Status (All == 3)",
                             "All DGA Status")
      }  else {
        title_name <- switch(L_n,
                             `1`  = "DGA Status 1",
                             `2`  = "DGA Status 2",
                             `3`  = "DGA Status 3",
                             `n1` = "DGA Status 2 or 3")
      }
      # input data order does not match plot, reordering:
      sampled <- sampled[,c(2,1,3)]
      default <- default[c(2,1,3)]
      TernaryPlot(point=DTn$labels[[1]],
                  atip=DTn$labels[[2]], btip=DTn$labels[[3]], ctip=DTn$labels[[4]],
                  alab=DTn$labels[[5]], blab=DTn$labels[[6]], clab=DTn$labels[[7]])
      for(idx in c(1:length(DTn$zones))) {
        TernaryPolygon(DTn$zones[[idx]], col=DTn$colours[[idx]], border='grey')
      }
      title(main = paste0("Duval Triangle ", DT_n, ": ", title_name))
      legend('topright',
             legend=unique(unlist(lapply(strsplit(names(DTn$zones), "_"), function(x) x[[1]]))),
             cex=0.8, bty='n', pch=21, pt.cex=1.8,
             pt.bg=unique(unlist(DTn$colours)))
      legend('topleft',
             legend=c("sampled value", "default value"),
             cex=0.8, bty='n', pch=c(".","*"), pt.cex=1.8)
      TernaryPoints(sampled, col=rgb(0,0,0,0.5), pch = ".") # Add points
      TernaryDensityContour(sampled, resolution = 50L) # Contour by point density
      TernaryPoints(default, col = "white", pch = "*") # Add the default location
      # Plots frequency of given diagnosis
      DTn <- paste0("DT", DT_n, "_L", L_n)
      plot(DI[names(DI) == DTn][[1]],
           main = paste0("Duval Triangle ", DT_n, ": ", title_name),
           xlab = "Diagnosis",
           ylab = "Proportion of Samples with Diagnosis")
    }
  }

  if((PP$diagnosis == 1) | (PP$diagnosis > 2)) { # RR | IEC 1|2
    Summarise_Abs_n <- function(Abs_n, colour_scale, x_scale, y_scale, df, df_0, L1_mask, L2_mask, L3_mask, skip) {
      if(Abs_n == "RR") {
        Abs_n_zones <- Define_RR()
        Abs_n_name <- "Rogers Ratios Diagnosis"
      }
      if(Abs_n == "IEC1") {
        Abs_n_zones <- Define_IECn("1")
        Abs_n_name <- "IEC Table 1 Diagnosis"
      }
      if(Abs_n == "IEC2") {
        Abs_n_zones <- Define_IECn("2")
        Abs_n_name <- "IEC Table 2 Diagnosis"
      }
      Plot_Abs_n(Abs_n_zones, Abs_n, Abs_n_name, colour_scale, x_scale, y_scale, df, df_0, "0", skip)
      if(skip == 0) {
        if(sum(L1_mask)>0) Plot_Abs_n(Abs_n_zones, Abs_n, Abs_n_name, colour_scale, x_scale, y_scale, df[rep(L1_mask,2),], df_0, "1", 0)
        # if(sum(L2_mask)>0) Plot_Abs_n(Abs_n_zones, Abs_n, Abs_n_name, colour_scale, x_scale, y_scale, df[rep(L2_mask,2),], df_0, "2", 0)
        # if(sum(L3_mask)>0) Plot_Abs_n(Abs_n_zones, Abs_n, Abs_n_name, colour_scale, x_scale, y_scale, df[rep(L3_mask,2),], df_0, "3", 0)
        if(sum(!L1_mask)>0) Plot_Abs_n(Abs_n_zones, Abs_n, Abs_n_name, colour_scale, x_scale, y_scale, df[rep(!L1_mask,2),], df_0, "n1", 0)
      }
    }

    Plot_Abs_n <- function(Abs_D, Abs_n, Abs_n_name, colour_scale, x_scale, y_scale, sampled, default, L_n, skip) {
      if(L_n == "0") {
        title_name <- "All DGA Status"
        if(skip > 0) title_name <- paste0(title_name, " (All == ", skip, ")")
      }  else {
        title_name <- switch(L_n,
                             `1`  = "DGA Status 1",
                             `2`  = "DGA Status 2",
                             `3`  = "DGA Status 3",
                             `n1` = "DGA Status 2 or 3")
      }
      suppressWarnings(suppressMessages(print(ggplot(sampled) +
          geom_rect(data = Abs_D, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=Abs_D$label),alpha=0.3) +
          geom_point(aes(x = x, y = y, colour = D)) +
          geom_point(data = default, aes(x = x, y = y), fill = 'white', colour = 'white', shape=8) +
          facet_grid(~ name) + theme(legend.position="bottom") +
          scale_x_log10(breaks = x_scale$breaks, limits = x_scale$limits) +
          scale_y_log10(breaks = y_scale$breaks, limits = y_scale$limits) +
          scale_colour_manual(values = colour_scale, drop = FALSE) +
          scale_fill_manual(values = colour_scale, drop = FALSE) +
          guides(fill=guide_legend(title="Diagnosis", override.aes=list(colour=colour_scale)), colour="none") +
          ggtitle(paste0(Abs_n_name, ": ", title_name), subtitle = "* = sampled value, • = default value") +
          xlab("Gas Ratio") + ylab("C2H2/C2H4"))))
      # Plots frequency of given diagnosis
      Abs_n <- paste0(Abs_n, "_L", L_n)
      plot(DI[names(DI) == Abs_n][[1]],
           main = paste0(Abs_n_name, ": ", title_name),
           xlab = "Diagnosis",
           ylab = "Proportion of Samples with Diagnosis")
    }
  }

  if((PP$diagnosis == 1) | (PP$diagnosis == 2)) { # DT
    Summarise_DTn("1", DI$DT1_0_data, DI$DT1_default_data,
                  SC$L1_mask, SC$L2_mask, SC$L3_mask, skip)
    Summarise_DTn("4", DI$DT4_0_data, DI$DT4_default_data,
                  SC$L1_mask, SC$L2_mask, SC$L3_mask, skip)
    Summarise_DTn("5", DI$DT5_0_data, DI$DT5_default_data,
                  SC$L1_mask, SC$L2_mask, SC$L3_mask, skip)
  }

  if((PP$diagnosis == 1) | (PP$diagnosis > 2)) { # RR|IEC: Abs. Ratios [C2H2:C2H4, CH4:H2, C2H4:C2H6]
    # These are the default values
    df_0 <- data.frame(y = unlist(rep(DI$Abs_default_data[1],2)),
                       x = unlist(c(DI$Abs_default_data[2],
                                    DI$Abs_default_data[3])),
                       name = c("CH4/H2", "C2H4/C2H6"))
    # These are the sampled points
    df <- data.frame(y = rep(DI$Abs_0_data[1,],2),
                     x = c(DI$Abs_0_data[2,],
                           DI$Abs_0_data[3,]),
                     name = c(rep("CH4/H2", length(DI$Abs_0_data[1,])),
                              rep("C2H4/C2H6", length(DI$Abs_0_data[1,]))))
    # Just for visuals - moving points out of bounds back into bounds.
    df_0$x[df_0$x > 10.5] <- 10.5
    df_0$y[df_0$y > 10.5] <- 10.5
    df_0$x[df_0$x < 0.01] <- 0.005
    df_0$y[df_0$y < 0.01] <- 0.005
    df$x[df$x > 10.5] <- 10.5
    df$y[df$y > 10.5] <- 10.5
    df$x[df$x < 0.01] <- 0.005
    df$y[df$y < 0.01] <- 0.005
  }

  if((PP$diagnosis == 1) | (PP$diagnosis == 3)) { # RR
    colour_scale <- c('darkgrey','green','yellow','mediumorchid',
                      'tomato','red','darkred','black')
    x_scale <- list(breaks = c(0.10, 0.5, 1.0, 3.0, 10), limits = c(0.005, 11))
    y_scale <- list(breaks = c(0.10, 0.5, 1.0, 3.0, 10), limits = c(0.005, 11))
    names(colour_scale) <- c('.Unknown','C0: Normal','C1: PD','C2: Arcing',
                             'C3: Thermal Low','C4: Thermal<700C','C5: Thermal>700C','Truncated')
    df_0$D <- unlist(rep(DI$RR_D_default[1],2))
    levels(df_0$D) <- names(colour_scale) # Make labels match
    df$D <- unlist(rep(DI$RR_D,2))
    levels(df$D) <- names(colour_scale) # Make labels match
    # Plots relevant plots
    Summarise_Abs_n("RR", colour_scale, x_scale, y_scale, df, df_0,
                    SC$L1_mask, SC$L2_mask, SC$L3_mask, skip)
  }

  if((PP$diagnosis == 1) | (PP$diagnosis == 4)) { # IEC 1|2
    # IEC 1
    colour_scale <- c('darkgrey','yellow','lightpink','darkorchid',
                      'tomato','red','darkred','black')
    x_scale <- list(breaks = c(0.10, 0.2, 0.5, 1.0, 2.0, 4.0, 10), limits = c(0.005, 11))
    y_scale <- list(breaks = c(0.01, 0.1, 0.2, 0.6, 1.0, 2.5, 10), limits = c(0.005, 11))
    names(colour_scale) <- c('.Unknown','1.PD','2.Arcing Low','3.Arcing High',
                             '4.Thermal<300C','5.Thermal<700C','6.Thermal>700C','Truncated')
    df_0$D <- unlist(rep(DI$IEC1_D_default[1],2))
    levels(df_0$D) <- names(colour_scale) # Make labels match
    df$D <- unlist(rep(DI$IEC1_D,2))
    levels(df$D) <- names(colour_scale) # Make labels match
    # Plots relevant plots
    Summarise_Abs_n("IEC1", colour_scale, x_scale, y_scale, df, df_0,
                    SC$L1_mask, SC$L2_mask, SC$L3_mask, skip)
    # IEC 2
    colour_scale <- c('darkgrey','yellow','mediumorchid','red','black')
    names(colour_scale) <- c('.NA','1.PD','2.Arcing','3.Thermal','Truncated')
    # Only need "CH4/H2", Keeping both here but only plotting needed one.
    mask <- df$name == "CH4/H2"
    df_0$D <- unlist(rep(DI$IEC2_D_default[1],2))
    levels(df_0$D) <- names(colour_scale)
    df$D <- unlist(rep(DI$IEC2_D,2))
    levels(df$D) <- names(colour_scale)
    # Plots relevant plots
    Summarise_Abs_n("IEC2", colour_scale, x_scale, y_scale, df[mask,], df_0[1,],
                    SC$L1_mask, SC$L2_mask, SC$L3_mask, skip)
  }
}
