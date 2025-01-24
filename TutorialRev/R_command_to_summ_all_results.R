source("./custom_R_functions/sum_and_plot_traderpros_result.R")



results <- sum_and_plot_traderpros_result(
     out_dir = "../outputfIGURESpy2",
     in_dir = "../outputsTESTpy2",
     prefix_used_in_traderpros = "CicTraderHiSSE",
     burn = 0.1,
     out_images_format = "both",
     path_to_custom_functions = "./custom_R_functions"
   )

results
