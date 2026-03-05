plot_gg_tree_alex <- function (
    tree,
    common_name = NA,
    add_x = TRUE,
    add_tip_label = "short_tip_lab",
    add_bootstrap = FALSE,
    add_branch_length = FALSE,
    add_branch_id = FALSE,
    branch_text_param = NULL,
    branch_text_size = 3,
    only_shared_branches = FALSE,
    plot_margin = grid::unit(c(10, 100, 10, 10), "points"),
    branch_color_param = NULL,
    branch_colors = NA,
    add_title = TRUE,
    legend_pos = NULL
) {
  # --- Safe common_name handling
  if (is.na(common_name)) {
    if (exists("find_common_sample_string", mode = "function")) {
      common_name <- find_common_sample_string(tree)
    } else {
      common_name <- NULL
    }
  }
  
  # --- Tip labels (needs stringr)
  if (!is.null(common_name)) {
    tree@data$short_tip_lab <- stringr::str_remove(tree@data$tip.label, common_name)
  } else {
    tree@data$short_tip_lab <- tree@data$tip.label
  }
  
  # --- Only shared branches
  if (only_shared_branches) {
    tree@data[tree@data$n_samples == 1, setdiff(colnames(tree@data), "node")] <- NA
    add_tip_label <- NA
    plot_margin <- grid::unit(c(10, 10, 10, 10), "points")
  }
  
  # --- Colors fallbacks (in case your palettes aren't defined)
  if (!exists("dist_cols", inherits = TRUE))    dist_cols    <- scales::hue_pal()(22)
  if (!exists("dist_cols31", inherits = TRUE))  dist_cols31  <- scales::hue_pal()(31)
  if (!exists("grad_cols", inherits = TRUE))    grad_cols    <- c("white", "black")
  
  # --- Base tree
  if (!is.null(branch_color_param)) {
    gg <- ggtree::ggtree(
      tree,
      branch.length = "branch_length",
      mapping = ggplot2::aes(color = rlang::sym(branch_color_param))
    ) + ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(plot.margin = plot_margin)
    
    if (length(branch_colors) == 1 && is.na(branch_colors)) {
      # auto colors
      if (any(c("factor", "character") %in% class(dplyr::pull(tree@data, !!branch_color_param)))) {
        classes <- length(unique(tree@data[, branch_color_param]))
        if (classes <= 22)      cols <- dist_cols
        else if (classes <= 31) cols <- dist_cols31
        else if (classes <= 50) cols <- dist_cols31
        else                    cols <- (scales::hue_pal())(classes)
        gg <- gg + ggplot2::scale_color_manual(values = cols)
      } else {
        gg <- gg + ggplot2::scale_color_gradientn(colors = grad_cols)
      }
    } else {
      # user colors provided
      if (any(c("factor", "character") %in% class(tree@data[, branch_color_param]))) {
        gg <- gg + ggplot2::scale_color_gradientn(colors = branch_colors)
      } else {
        gg <- gg + ggplot2::scale_color_manual(values = branch_colors)
      }
    }
  } else {
    gg <- ggtree::ggtree(tree, branch.length = "branch_length") +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(plot.margin = plot_margin)
  }
  
  # --- Optional annotations
  if (add_branch_id) {
    gg <- gg + ggplot2::geom_text(ggplot2::aes(x = branch, label = branch_id),
                                  nudge_y = 0.5, size = 3, na.rm = TRUE)
  }
  if (add_branch_length) {
    gg <- gg + ggplot2::geom_text(ggplot2::aes(x = branch, label = branch_length),
                                  nudge_y = -0.5, size = 2, na.rm = TRUE)
  }
  
  # --- Bootstrap as burgundy dots on the node
  if (add_bootstrap) {
    gg <- gg +
      ggtree::geom_point2(
        aes(subset = !is.na(n_boot), fill = as.numeric(n_boot)), # force numeric
        shape = 21, color = "black", size = 3, stroke = 0.3, na.rm = TRUE
      ) +
      ggplot2::scale_fill_gradient(
        low = "white", high = "#800020", limits = c(0, 100),
        name = "Bootstrap"
      )
  }
  
  # --- Tip labels
  if (!is.na(add_tip_label)) {
    if (any(class(add_tip_label) == "character")) {
      gg <- gg + ggtree::geom_tiplab(ggplot2::aes(label = rlang::sym(add_tip_label)),
                                     size = 3, color = "grey70")
    } else if (is.logical(add_tip_label) && add_tip_label) {
      gg <- gg + ggtree::geom_tiplab(size = 3, color = "grey70")
    }
  }
  
  # --- Optional branch text
  if (!is.null(branch_text_param)) {
    gg <- gg + ggplot2::geom_text(
      ggplot2::aes(x = branch, label = rlang::sym(branch_text_param)),
      nudge_y = 0.1, size = branch_text_size, color = "black"
    )
  }
  
  # --- X axis (safe fallback if helper is missing)
  if (isTRUE(add_x)) {
    if (exists("add_x_axis", mode = "function")) {
      gg <- gg + add_x_axis()
    } else {
      # minimal x-axis so the plot keeps its scale without error
      gg <- gg + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.05)))
    }
  }
  
  # --- Title
  if (is.logical(add_title) && add_title) {
    gg <- gg + ggplot2::labs(title = common_name)
  } else if (is.character(add_title)) {
    gg <- gg + ggplot2::labs(title = add_title)
  }
  
  # --- Legend position
  if (!is.null(legend_pos)) {
    gg <- gg + ggplot2::theme(legend.position = legend_pos)
  }
  
  gg
}
