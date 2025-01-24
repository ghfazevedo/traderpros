#' Summary Statistics for Trace Variables
#'
#' This function calculates summary statistics for a specified set of variables within an MCMC trace file or object.
#' Either `trace_file` or `trace_object` must be specified.
#'
#' @param trace_file Path to the trace file. Either `trace_file` or `trace_object` must be specified.
#' @param trace_object An existing trace object created with `readTrace()`. Either `trace_file` or `trace_object` must be specified.
#' @param variables A character vector of variable names to summarize. Must be present in the trace file or trace object.
#' @param variable_match A vector of string patterns to match variable names in the trace data using `grep`, used if `variables` is not specified.
#' @param out_name (Optional) The name of the output file to save the summary as a CSV.
#' @param burn_in The burn in percentage if reading trace file. Default is 0.1
#'
#' @return A data frame containing summary statistics of the specified trace variables.
#' @export
#' 

summary_stats_from_trace <- function(trace_file = NULL,
                                    trace_object = NULL,
                                    variables = NULL,
                                    variable_match = NULL,
                                    out_name = NULL,
                                    burn_in = 0.1) {
  
  # Load trace data from file if not provided as an object
  if (!is.null(trace_file)) {
    traceModel <- readTrace(trace_file, burnin = burn_in)
  } else if (!is.null(trace_object)) {
    traceModel <- trace_object
  } else {
    stop("Either 'trace_file' or 'trace_object' must be specified.")
  }
  
  # Determine variables to summarize
  if (!is.null(variable_match)) {
    # Apply grep to each pattern in variable_match, unlist, and take unique values
    matched_vars <- unique(unlist(lapply(variable_match, function(pattern) {
      grep(pattern, names(traceModel[[1]]), value = TRUE)
    })))
    
    # Check if any variables were matched
    if (length(matched_vars) == 0) {
      stop("No variables matched the specified patterns.")
    }
    variables <- matched_vars
  } else if (is.null(variables)) {
    stop("Either 'variables' or 'variable_match' must be specified.")
  }
  
  # Summarize the specified trace variables
  RatesSummary <- as.data.frame(summarizeTrace(trace = traceModel, vars = variables))
  
  # Set column names
  names(RatesSummary) <- variables
  
  # Optionally save output
  if (!is.null(out_name)) {
    write.csv(as.data.frame(t(RatesSummary)), out_name, row.names = TRUE)
  }
  
  # Return the summary data frame
  return(RatesSummary)
}
