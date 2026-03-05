calculate_product_sensitivity <- function(samples, df_sensitivity) {
  # Split the 'samples' string by '|'
  sample_list <- strsplit(samples, "\\|")[[1]]
  
  # Retrieve sensitivities for these samples
  sensitivities <- df_sensitivity$sensitivity[df_sensitivity$Sample_name %in% sample_list]
  
  # Calculate the product of (1 - sensitivity)
  1 - prod(1 - sensitivities)
}