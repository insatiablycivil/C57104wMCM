#' Screening_Plot: Plotting outputs from Screening function for validation.
#'
#' This is internal function to plot generated screening outputs for validation.
#'
#' @param SC List: Outputs from Screening function.
#' @param PP List: Outputs from Preprocess function.
#' @examples
#' # Not recreating required inputs to run this example
#' # Screening_Plot(SC = screened,
#' #                PP = preprocessed)
#' @export
Screening_Plot <- function(SC, PP) {
  # For viewing the data for validation purposes.
  # SC has screening output
  # PP (preprocessed) has various parameters
  # Looking at:
  #   1: gas table inputs and screening outputs
  #   2: gas screening outputs and combined screening output

  Upper_Plot_Cont <- function(data, mapping, ...) {
    # TODO: This seems dumb, should be a better way
    x <- unlist(lapply(mapping, all.vars))
    i_1 <- match(x[1], names(data))
    i_2 <- match(x[2], names(data))
    w_1 <- data[,i_1]
    w_2 <- data[,i_2]
    if((length(w_1)<3) | (length(w_2)<3)) {
      p <- ggally_text("Too Few Samples!")
    } else {
      test <- stats::cor.test(w_1, w_2)
      cor_text <- paste0("Cor=", round(as.numeric(test$estimate), 4))
      if(test$p.value <= 0.05) cor_text  <- paste0(cor_text,"*")
      if(test$p.value <= 0.01) cor_text  <- paste0(cor_text,"*")
      if(test$p.value <= 0.001) cor_text <- paste0(cor_text,"*")
      p <- ggally_text(cor_text)
    }
    p
  }

  Diag_Plot_Cont <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_histogram(alpha = 0.5, bins = bin_count)
    p
  }

  Diag_Plot_Discrete <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_bar() + scale_x_discrete(drop = FALSE)
    p
  }

  Diag_Plot_Discrete_Alt <- function(data, mapping, ...) {
    # TODO: This seems dumb, should be a better way
    x <- unlist(lapply(mapping, all.vars))
    i <- match(x, names(data))
    w <- data[,i]
    p <- ggplot(data = data, mapping = mapping) +
      geom_bar(aes(fill = w, colour = w), stat="count") + scale_x_discrete(drop=FALSE)
    p
  }

  Lower_Plot_Cont <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(alpha = 0.2) +
      geom_smooth(method = "lm", color = "black", ...)
    p
  }

  Lower_Plot_Discrete <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_col() + scale_x_discrete(drop=FALSE)
    p
  }

  Lower_Plot_Combo <- function(data, mapping, ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_jitter(alpha = 0.2) +
      scale_y_discrete(drop = FALSE)
    p
  }

  Lower_Plot_Discrete_Alt <- function(data, mapping, ...) {
    # TODO: This seems dumb, should be a better way
    x <- unlist(lapply(mapping, all.vars))
    i_1 <- match(x[1], names(data))
    i_2 <- match(x[2], names(data))
    w_1 <- data[,i_1]
    w_2 <- data[,i_2]
    data = data.frame(x = w_1, y = w_2)
    p <- ggplot(data = data, aes(x=w_1)) +
      geom_bar(aes(fill = w_2, colour = w_2), position="stack", stat="count") +
      scale_x_discrete(drop = FALSE)
    p
  }

  Lock_Fill_Scale <- function(...){
    ggplot2:::manual_scale(
      'fill',
      values = setNames(c('green', 'orange', 'red', 'black'), c("1","2","3","black")),
      drop = FALSE,
      ...
    )
  }

  Lock_Colour_Scale <- function(...){
    ggplot2:::manual_scale(
      'colour',
      values = setNames(c('green', 'orange', 'red', 'black'), c("1","2","3","black")),
      drop = FALSE,
      ...
    )
  }

  d <- dim(SC$T12_vals)[1]
  g <- dim(SC$T12_vals)[2]
  n <- dim(SC$T12_vals)[3]

  # Get Table Metrics, Per Gas Status (combined to 1:3), and Combined Status (combined to 1:3)
  screening_vals <- abind(SC$T12_vals, SC$T3_vals, SC$T4_vals,
                          SC$GL1_mask+(2*SC$GL2_mask)+(3*SC$GL3_mask),
                          t(array(rep(SC$L1_mask+(2*SC$L2_mask)+(3*SC$L3_mask),g), dim=c(n, g))),
                          along = 1)
  # Plotting 1: Gas Table Inputs and Status Outputs
  bin_count <- ceiling(log(n)*2)
  for(idx in c(1:g)) {
    df <- as.data.frame(t(screening_vals[,idx,]))
    colnames(df) <- c("Abs Value (T1|2)", "Delta (T3)", "Gassing Rate (T4)", "Gas Status", "Combined Status")
    df$`Gas Status` <- factor(df$`Gas Status`, levels = c(1,2,3))
    df$`Combined Status` <- factor(df$`Combined Status`, levels = c(1,2,3))
    title_name <- paste0("Correlations Across Screening Metrics and Outputs: Gas ", PP$col_names[idx])
    suppressWarnings(suppressMessages(print(ggpairs(df,
            mapping = ggplot2::aes(color = `Gas Status`, fill = `Gas Status`),
            upper = list(continuous = wrap(Upper_Plot_Cont)),
            diag = list(continuous = wrap(Diag_Plot_Cont),
                        discrete = wrap(Diag_Plot_Discrete)),
            lower = list(continuous = wrap(Lower_Plot_Cont),
                         discrete = wrap(Lower_Plot_Discrete),
                         combo = wrap(Lower_Plot_Combo)),
            title = title_name) + Lock_Fill_Scale() + Lock_Colour_Scale())))
  }

  # Get Per Gas Status (combined to 1:3) and Combined Status (combined to 1:3)
  if(g > 1) {
    df <- as.data.frame(t(screening_vals[4,,]))
    df$V8 <- SC$L1_mask+(2*SC$L2_mask)+(3*SC$L3_mask)
    colnames(df) <- c(PP$col_names, "Combined Status")
    df[] <- lapply(df, function(x) factor(as.character(x), levels = c(1,2,3)))
    # Plotting 2: Gas Status Outputs and Combined Status
    suppressWarnings(suppressMessages(print(ggpairs(df,
                  upper = "blank",
                  diag = list(discrete = wrap(Diag_Plot_Discrete_Alt)),
                  lower = list(discrete = wrap(Lower_Plot_Discrete_Alt)),
                  title = "Correlations Across Screening Outputs") +
            Lock_Fill_Scale() + Lock_Colour_Scale() +
            ylim(0, n))))
  }
}
