# patient colors from the HSCT paper. 7 sets of 1 dark, 1 light color
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

# gradient of 5 colors
grad_cols <- c(gplots::colorpanel(n = 50, low = '#2f2bad', mid = '#93cfe2', high = '#fffaaa'),
             gplots::colorpanel(n = 50, low = '#fffaaa', mid = "#ffad66", high = '#d60000'))

# colors I use for aging/treatment signatures
signature_cols = c("SBS18" = "#86b7bd",
                   "SBS1" = "grey50", 
                   "SBS5" = "grey70",
                   "HSPC" = "grey90",
                   "SBS87" = "#B6244F",
                   "SBS31" = "#f79d5c",
                   "SBS35" = "#bd692d",
                   "SBS17a" = "#e8ed7e",
                   "SBS17b" = "#9fa33c",
                   "SBSA" = "#a269c2",
                   "SBSB" = "#6bb388",
                   "SBSC" = "#4d5d91",
                   'unexplained' = '#e8d0ac')

# presets for R panels
dark_cols = palette('Dark2')

paired_cols = palette('Paired')

accent_cols = palette('Accent')

if (FALSE) {
  '
  A few people have been asking me about color schemes recently. Id thought to share the tools I use.
  
  First, the R color palettes are quite nice, certainly for a few samples, 
  or highly distinguishable colors. From R.4 onwards you can use "palette()" to generate different ones, 
  e.g. "palette("Dark2")" that Sjors also uses. 
  https://blog.r-project.org/2019/11/21/a-new-palette-for-r/
  
  To create gradients I use gplots::color_panel, with a low, mid and high color 
  (and combining two of these to get a gradient of 5 colors)
  
  To generate small color schemes: https://coolors.co/ is really nice. 
  You can play around with colors, and if you pick a few, 
  it will generate suiting colors for you. 
  Unfortunately, only 5 are free now, but if you want more, just lock 4 and generate a few.
  
  If you want to use the standard colors for journals like Nature or Lancet, use the "ggsci" package 
  (e.g. ggsci::pal_npg()(10) to get the 10 Nature colors).
  '
}
  
