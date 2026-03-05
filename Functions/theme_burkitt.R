# plotting theme for CHemALL sequencing data

# usage:
#source("~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R")
#theme_set(theme_CHemALL())

sixcols <- c('#BDD9BF',
             '#7785AC',
             '#D57A66',
             '#022F40',
             '#6B2737',
             '#0D5D56')

sixcols3 <- c('#B1D2B3','#BC7A8F','#8290BB','#54BFB7','#f3d08c','#D57A66')

extra_cols <- c('#0B5563','#472836')
extra_cols2 <- c('#dc9c38', '#3E5496')

### merge for extended palettes
sevencols <- c('#474647', sixcols3)
eightcols <- c(sixcols3, extra_cols2)
eightcols_withref <- c('#474647', eightcols)

chemall_colors <- eightcols_withref
names(chemall_colors) <- c("healthy \nreference", "PB13467",
                           "PB17708","PB22514","PB14041","PB20885", "PB08163" , "PB06238","PB11226"  )

### colors for QC metric and signatures
paired_six <- c('#54BFB7','#0A9086','#BC7A8F','#8E2043','#8290BB','#3E5496')

paired_six2 <- c('#f3d08c','#dc9c38','#BC7A8F','#8E2043','#8290BB','#3E5496') # alternative reds e27d8d 8d2f53

sig_cols <- c('beige','#0A9086','#b3b3b3', '#e5e5e5','#A62639', paired_six2[5:6], '#d55e00','#54BFB7', '#d55e00','#54BFB7',paired_six2[1:2])
names(sig_cols) <- c('unexplained','SBS18','SBS5','SBS1','HSPC','PTA_Middelkamp','PTA_Luquette','ResolveOMEv1','ResolveOMEv2', 'gtPTAv1','gtPTAv2', 'SBS19','SBS87')

### colors for mutation types
mut_colors <- c(
  "#2EBAED", "#000000", 
  "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE", "white"
)
names(mut_colors) <- c('C>A','C>G','C>T','T>A','T>C','T>G', 'rest')

## optimized version of eight
palette_8_cht <- c('#3E5496', '#F3D08C', '#A56CC1', '#54BFB7', '#A88F2A', '#8290BB', '#D57A66', '#B1D2B3')

##### historic Boxtel lab colors (use for healthy data if needed)

hsct_cols = c("#bc3c29", '#f0a9a8',"#e08f3e",'#fcd1ae','#f7e454', "#faf6d4", "#177354", "#76b874", "#45ad9c", '#8ae6e6', "#593C8F", "#CAA8F5", "#700657", "#AF4D98")

# distinguishable colors, from aesthetically nice to less nice, but with increasingly more colors
dist_cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
               '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
               '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
               '#000075', '#808080', '#b8ccd1', '#000000')
dist_cols31 <- c('#ff0000', '#7f2200', '#ffe1bf', '#474d00', '#ccff00', '#00593c', '#00ccff', '#002e73', '#220033',
                 '#ff00aa', '#ff4059', '#330000', '#ffa280', '#8c5e00', '#fbffbf', '#71b359', '#00ffcc', '#005580',
                 '#0044ff', '#cc00ff', '#7f0044', '#806460', '#ff8800', '#ffcc00', '#64664d', '#00ff22', '#7ca69d',
                 '#0088ff', '#d0bfff', '#9553a6', '#ffbfd9')
dist_cols50 <- c('#808080', '#c0c0c0', '#2f4f4f', '#556b2f', '#8b4513', '#2e8b57', '#7f0000', '#191970',
                 '#006400', '#808000', '#5f9ea0', '#b8860b', '#4682b4', '#d2691e', '#9acd32', '#cd5c5c',
                 '#00008b', '#32cd32', '#8fbc8f', '#8b008b', '#d2b48c', '#ff4500', '#ffa500', '#ffd700',
                 '#0000cd', '#00ff00', '#9400d3', '#00ff7f', '#4169e1', '#dc143c', '#00ffff', '#00bfff',
                 '#adff2f', '#ff6347', '#ff00ff', '#db7093', '#f0e68c', '#ffff54', '#6495ed', '#dda0dd',
                 '#87ceeb', '#ff1493', '#ffa07a', '#ee82ee', '#98fb98', '#7fffd4', '#ff69b4', '#fffacd', '#e0ffff', '#ffb6c1')






theme_CHemALL <- function (base_size = 12, base_family = "")
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(text =  element_text(size =  10, color = 'black'),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.3, color = "grey80"),
          panel.background = element_blank(),
          axis.line	= element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          legend.key = element_rect(fill = NA, color = NA)
    ) +
    theme(text =  element_text(size =  7, color = 'black'),
          axis.text = element_text(size = 7, colour = "black")) +
    theme(legend.key.size = unit(0.6, "lines"),              
          legend.text = element_text(margin = margin(l = 4)),   
          legend.text.align = 0, 
          legend.spacing.y = unit(0.2, "lines"),    
          legend.spacing.x = unit(0.2, "lines") )
  
}

ggEmptyAxis <- function(xy = 'x', line = F, text = T, ticks = T, title = T, preset = 'allbutline') {
  thm = theme()
  if (grepl('x', xy)) {
    if (line) { thm = thm + theme(axis.line.x = element_blank()) }
    if (text) { thm = thm + theme(axis.text.x = element_blank()) }
    if (ticks) { thm = thm + theme(axis.ticks.x = element_blank()) }
    if (title) { thm = thm + theme(axis.title.x = element_blank()) }
  }
  if (grepl('y', xy)) {
    if (line) { thm = thm + theme(axis.line.y = element_blank()) }
    if (text) { thm = thm + theme(axis.text.y = element_blank()) }
    if (ticks) { thm = thm + theme(axis.ticks.y = element_blank()) }
    if (title) { thm = thm + theme(axis.title.y = element_blank()) }
  }
  thm
}

ggTextAxisRotate <- function(xy = 'x', angle = 45) {
  thm = theme()
  if (grepl('x', xy)) {
    if (angle == 45) { thm = thm + theme(axis.text.x = element_text(angle = angle, hjust = 1)) }
    if (angle == 90) { thm = thm + theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 0.5)) }
  }
  if (grepl('y', xy)) {
    if (angle == 45) { thm = thm + theme(axis.text.y = element_text(angle = angle, hjust = 1)) }
    if (angle == 90) { thm = thm + theme(axis.text.y = element_text(angle = angle, hjust = 1, vjust = 0.5)) }
  }
  thm
}

facetnested_theme <- theme(
  strip.placement = 'outside',
  strip.background = element_rect(fill = "gray90", color = 'white', size = 0.5),
  strip.text.x = element_text(size = 5, margin = margin(t = 1, b = 1), color = "black")
)

profile96_theme <- theme(
  strip.text.x = element_text(size = 6),  # Adjust x-axis facet labels size
  strip.text.y = element_text(size = 6),    # smaller y facet labels
  plot.title = element_text(size = 8),   # Adjust y-axis facet labels size
  legend.title = element_text(size = 8),        # smaller legend title text
  legend.text = element_text(size = 6),         # smaller legend labels text
  legend.key.size = unit(0.5, "lines"),          # smaller legend boxes
  axis.title.x = element_text(size = 8),
  axis.title.y = element_text(size = 8),
  axis.text.x = element_text(size = 4),
  axis.text.y = element_text(size = 7)
)

collect_legend <- function(input_plot) {
  input_plot <- lapply(input_plot, function(p) p + theme(plot.margin = margin(5, 5, 5, 5)))
  output_plot <- wrap_plots(input_plot, guides = 'collect') & theme(legend.position = 'right')
  return(output_plot)
}

CHemALL_style <- list(
  theme_CHemALL(),
  ggTextAxisRotate(),
  scale_fill_manual(values = brewer.pal(n = 10, name = "Paired"))
)


