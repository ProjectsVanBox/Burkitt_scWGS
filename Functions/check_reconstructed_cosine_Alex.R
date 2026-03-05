check_reconstructed_cosine_Alex <- function (contribution, mut_matrix, signatures, tree) 
{
  recon = signatures %*% contribution
  recon_cosine = diag(cos_sim_matrix(recon, mut_matrix[, colnames(contribution)]))
  tibble(recon = recon_cosine, n_mut = tree@data$branch_length[match(colnames(contribution), 
                                                                     tree@data$branch_id)], branch = colnames(contribution)) %>% 
    ggplot(aes(x = log10(n_mut + 1), y = recon)) + geom_point() + 
    geom_hline(yintercept = 0.80, lty = 3) + geom_vline(xintercept = log10(201), 
                                                        lty = 3) + ggrepel::geom_text_repel(data = ~subset(.x, 
                                                                                                           n_mut > 200 & recon < 0.80), aes(label = branch)) + labs(x = "log10 # of mutations", 
                                                                                                                                                                    y = "cosine sim. reconstructed vs original profile")
}