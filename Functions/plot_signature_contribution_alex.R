plot_tree_contribution_bars_new_alex <- function(
    tree, signatures = NULL, contribution = NULL, mut_matrix = NULL,
    signature_colors = NULL, scaling = 1, remove_min = 20, title,
    bar_height = 0.03,
    branch_color = NULL,        # new optional arg
    x_name = "Mutation burden (SNVs)",
    x_limits = NULL,
    x_breaks = NULL,
    axis_text_x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1,
                                        margin = ggplot2::margin(t = 6)),
    axis_title_x = ggplot2::element_text(size = 12, vjust = 0.5,angle = 180, margin = ggplot2::margin(t = 6)),
    show_x_axis = TRUE,
    ...
) {
  if ((is.null(contribution) && is.null(signatures)) || (!is.null(contribution) && !is.null(signatures)))
    stop("please provide either a contribution matrix, or the signatures you want to fit on")
  if (!is.null(signatures) && is.null(mut_matrix))
    stop("please also provide a mutation matrix of the branches")
  
  if (!is.null(signatures)) {
    contribution <- fit_to_signatures_strict_tree(mut_matrix, signatures, remove_min = remove_min)
  }
  
  if (is.null(signature_colors)) {
    n_sig <- nrow(contribution)
    if (n_sig <= 22) cols <- dist_cols
    else if (n_sig <= 31) cols <- dist_cols31
    else if (n_sig <= 50) cols <- dist_cols31
    else cols <- (scales::hue_pal())(n_sig)
    signature_colors <- setNames(cols, rownames(contribution))
    signature_colors <- signature_colors[!is.na(names(signature_colors))]
  }
  if (is.null(names(signature_colors))) {
    names(signature_colors) <- rownames(contribution)
    signature_colors <- signature_colors[!is.na(names(signature_colors))]
  }
  
  fraction <- contribution %>%
    t() %>% as.data.frame() %>%
    tibble::rownames_to_column("node") %>%
    tidyr::pivot_longer(cols = -dplyr::starts_with("node"),
                        names_to = "sig", values_to = "contribution")
  
  fraction$node <- tree@data$node[match(fraction$node, tree@data$branch_id)]
  fractions_split <- split(fraction, fraction$node)
  
  bars <- lapply(fractions_split, function(df) {
    ggplot(df, aes(x = "", y = contribution, fill = sig)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = signature_colors) +
      theme_void() + theme(legend.position = "none")
  })
  
  nodes_slct <- fractions_split %>% sapply(function(x) dplyr::pull(x, node)[1]) %>% unname()
  plot <- plot_gg_tree(tree, add_tip_label = FALSE, add_title = title, ...)
  
  # --- Styling additions
  if (!is.null(branch_color)) {
    plot <- plot + ggtree::geom_tree(color = branch_color)
  }
  if (isTRUE(show_x_axis)) {
    plot <- plot +
      ggtree::theme_tree2() +
      ggplot2::scale_x_continuous(
        name   = x_name,
        limits = x_limits,
        breaks = x_breaks,
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::theme(
        axis.text.x  = axis_text_x,
        axis.title.x = axis_title_x
      )
  }
  
  max_x <- max(layer_scales(plot)$x$range$range)
  bar_lengths <- tree@data$branch_length / max_x
  bar_lengths <- bar_lengths[match(nodes_slct, tree@data$node)] * scaling
  
  for (i in seq_along(bars)) {
    plot <- ggtree::inset(
      tree_view = plot, insets = bars[i],
      width = bar_lengths[i], height = bar_height, x = "branch"
    )
  }
  
  fig_for_leg <- ggplot(fraction, aes(x = "", y = contribution, fill = sig)) +
    geom_col() + scale_fill_manual(values = signature_colors)
  legend <- cowplot::get_legend(fig_for_leg)
  
  cowplot::plot_grid(plot, legend, nrow = 1, rel_widths = c(1, 0.2))
}
