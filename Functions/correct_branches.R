correct_branches <- function(tree, callable_df){
  # get sensitivity df by taking fraction of total possible callable loci
  max_callable <- 2745186691
  callable_df$sensitivity <- callable_df$Callable_Loci / max_callable
  sensitivity_df <- callable_df[c('Sample_name','sensitivity')]
  
  print(sensitivity_df)
  
  # calculate the sensitivity per node/branch using the calculate_product_sensitivity function
  tree@data$product_sensitivity <- sapply(tree@data$samples, calculate_product_sensitivity, sensitivity_df)
  
  # correct the branch lengths
  tree@data$corr_branch_lengths <- tree@data$branch_length / tree@data$product_sensitivity
  
  print(tree@data$branch_length)
  print(tree@data$corr_branch_lengths)
  
  tree@data$corr_branch_lengths[is.na(tree@data$corr_branch_lengths)] <- 0 # for merged germline sample
  
  # store also in the phylo object, required for downstream analyses based on the phylo object
  tree@phylo$edge.length <- tree@data[order(tree@data$node),]$corr_branch_lengths
  
  # store in branch_length column for cellphyplotting use
  tree@data$branch_length <- tree@data$corr_branch_lengths
  
  print(tree@phylo$edge.length)
  
  return(tree)
}