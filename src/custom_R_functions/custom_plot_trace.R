#' Custom Plot Trace Function
#'
#' This function creates a customized density plot for specified variables in a trace file or trace object,
#' with flexible color and labeling options, and saves the output plot in PDF, PNG, or both formats.
#' It requires RevGadgets, ggplot2, coda.
#' Additionaly  RColorBrewer and viridis might be useful
#'
#' @param trace_file Path to the trace file. Either `trace_file` or `trace_object` must be specified.
#' @param trace_object An existing trace object, created with `readTrace()`. Either `trace_file` or `trace_object` must be specified.
#' @param variables A character vector of variable names to be plotted. Must be present in the trace file or trace object.
#' @param variable_match A string pattern to match variable names in the trace file using `grep`. Used if `variables` is not specified.
#' @param colors A character vector of colors for each variable, matching the length of `variables`.
#' @param labels_name A character vector for legend labels, matching the length of `variables`.
#' @param fig_out_name Path to save the output file, without extension.
#' @param out_format Output format, either "pdf", "png", or "both" for saving in both formats.
#' @param burn_in The burn in percentage if reading trace file. Default is 0.1 
#'
#' @details This function checks if `trace_file` or `trace_object` is specified and reads the trace data accordingly.
#' If `variable_match` is provided, it uses `grep` to select matching variables from the trace.
#' It then customizes the plot with specified colors, labels, and adds a dashed line marking the median for each variable.
#' 
#' @return A ggplot object of the customized plot, displayed and saved to file.
#' 
#' @examples
#' @examples
#' custom_plot_trace(trace_file = "example_trace.csv", variables = c("var1", "var2"), ...)
#'

library(ggplot2)
library(RevGadgets)
library(coda)
#library(RColorBrewer)
#library(viridis)


custom_plot_trace <- function(trace_file = NULL,
                              trace_object = NULL,
                              variables,
                              variable_match = NULL,
                              colors,
                              labels_name,
                              fig_out_name,
                              out_format,
                              burn_in = 0.1) {
  
  # Check for input validity
  if (is.null(trace_file) && is.null(trace_object)) {
    stop("Either trace_file or trace_object must be specified.")
  }
  
  if (!is.null(trace_file) && !is.null(trace_object)) {
    stop("Specify only one of trace_file or trace_object.")
  }
  
  if (is.null(variable_match) && missing(variables)) {
    stop("Either variables or variable_match must be specified.")
  }
  
  if (!is.null(variable_match)) {
    if (is.null(trace_object)) {
      trace_object <- readTrace(trace_file, burnin = burn_in)
    }
    variables <- grep(variable_match, names(trace_object[[1]]), value = TRUE)
  }
  
  if (length(variables) != length(colors)) {
    stop("Length of colors must match length of variables.")
  }
  
  if (length(variables) != length(labels_name)) {
    stop("Length of labels_name must match length of variables.")
  }
  
  # Load the trace model
  if (!is.null(trace_file)) {
    traceModel <- readTrace(trace_file)
  } else {
    traceModel <- trace_object
  }
  
  # Generate the plot
  plotSpCompl <- plotTrace(trace = traceModel, vars = variables)
  ggplot_object <- plotSpCompl[[1]]
  
  # Change fill colors in the ggplot object
  for (i in seq_along(variables)) {
    ggplot_object$layers[[i + 1]]$aes_params$fill <- colors[i]  # Fill colors for ribbons
  }
  
  # Calculate medians
  median_values <- sapply(variables, function(var) {
    median(ggplot_object$data$value[ggplot_object$data$Variable == var])
  })
  
  # Create a vector for color mappings
  variables_colors_vector <- setNames(colors, variables)
  
  # Customize the plot
  ggplot_custom <- ggplot_object +
    scale_color_manual(values = variables_colors_vector,
                       name = NULL,
                       labels = labels_name) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(color = "black")) +
    labs(title = NULL)
  
  # Add median lines
  for (i in seq_along(median_values)) {
    ggplot_custom <- ggplot_custom +
      geom_vline(xintercept = median_values[i], color = colors[i], linetype = "dashed", linewidth = 0.5)
  }
  
  # Plot the customized ggplot object
  print(ggplot_custom)
  
  # Save the plot
  if ("pdf" %in% out_format) {
    ggsave(filename = paste0(fig_out_name, ".pdf"), plot = ggplot_custom, device = "pdf")
  }
  if ("png" %in% out_format) {
    ggsave(filename = paste0(fig_out_name, ".png"), plot = ggplot_custom, device = "png")
  }
  
  return(ggplot_custom)
}





