#' Plot RJ State-Dependent Probability Pie Chart with ggplot2
#'
#' @description This function creates a pie chart for state-dependent probabilities from a trace file or trace object.
#' It saves the chart as a specified file format (PNG, PDF, or both).
#'
#' @param trace_file Path to the trace file. Either `trace_file` or `trace_object` must be specified.
#' @param trace_object An existing trace object, created with `readTrace()`. Either `trace_file` or `trace_object` must be specified.
#' @param variable A variable name to be plotted. Must be present in the trace file or trace object.
#' @param colors A character vector of colors of length 2. Default is `c("black", "grey")`.
#' @param labels_name A character vector of labels for the legend of length 2. Default is `c("No", "Yes")`.
#' @param fig_out_name Path to save the output file, without extension. Defaults to `variable` name if not provided.
#' @param out_format Output format. Accepts `"pdf"`, `"png"`, or `"both"` to save in both formats.
#' @param burn_in The burn in percentage if reading trace file. Default is 0.1 
#'
#' @return A ggplot2 pie chart of the state-dependent probabilities.
#' 
#' @export
plot_rj <- function(trace_file = NULL,
                    trace_object = NULL,
                    variable,
                    colors = c("black", "grey"),
                    labels_name = c("No", "Yes"),
                    fig_out_name = NULL,
                    out_format = c("pdf", "png", "both"),
                    burn_in = 0.1) {
  
  # Argument checks
  if (is.null(trace_file) && is.null(trace_object)) {
    stop("Either `trace_file` or `trace_object` must be specified.")
  }
  if (!is.null(trace_file) && !is.null(trace_object)) {
    stop("Please provide only one of `trace_file` or `trace_object`.")
  }
  
  # Load trace model
  if (!is.null(trace_file)) {
    traceModel <- readTrace(trace_file, burnin = burn_in)
  } else {
    traceModel <- trace_object
  }
  
  # Check for valid variable
  if (!(variable %in% names(traceModel[[1]]))) {
    stop("The specified variable is not found in the trace data.")
  }
  
  # Set default fig_out_name if not provided
  if (is.null(fig_out_name)) {
    fig_out_name <- variable
  }
  
  # Validate colors and labels length
  if (length(colors) != 2) {
    stop("Colors must be a vector of length 2.")
  }
  if (length(labels_name) != 2) {
    stop("`labels_name` must be a vector of length 2.")
  }
  
  # Validate output format
  if (!out_format %in% c("pdf", "png", "both")) {
    stop("`out_format` must be one of 'pdf', 'png', or 'both'.")
  }
  
  
  # Generate probability data and label states
  # Count occurrences of 0 and 1, ensuring both are present
  counts <- table(factor(traceModel[[1]][[variable]], levels = c(0, 1)))
  
  # Compute percentages
  percentages <- round(100 * counts / sum(counts), 1)
  
  # Prepare data for ggplot
  data <- data.frame(
    State = labels_name, 
    Count = as.numeric(counts),
    Percentage = percentages
  )
  
  # Filter out missing categories
  data <- data[data$Count > 0, ]
  
  # Create the ggplot pie chart
  ggplot_pie <- ggplot(data, aes(x = "", y = Count, fill = State)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    scale_fill_manual(values = colors) +
    theme_void() +
    labs(title = variable, fill = "State") + # Add title and legend title
    geom_text(aes(label = paste0(State, ": ", Percentage.Freq, "%")),
              position = position_stack(vjust = 0.5), color = "white") +
    theme(legend.position = "none") # Optionally remove legend
  
  
  # Display the plot
  print(ggplot_pie)
  
  # Save plot as specified formats
  if (out_format == "png" || out_format == "both") {
    ggsave(filename = paste0(fig_out_name, ".png"), plot = ggplot_pie, width = 6, height = 6, dpi = 300)
  }
  
  if (out_format == "pdf" || out_format == "both") {
    ggsave(filename = paste0(fig_out_name, ".pdf"), plot = ggplot_pie, width = 6, height = 6)
  }
  
  # Return the ggplot object
  return(ggplot_pie)
}
