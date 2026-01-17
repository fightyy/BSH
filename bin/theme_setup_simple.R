###private theme
lwd_pt <- .pt * 72.27 / 96
base_size <- 7


my_theme <- function() {
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.5 / lwd_pt),
    axis.text.x = element_text(size = base_size, color = "black", 
                               lineheight = 0.9, 
                               margin = margin(t = 0.05, r = 0, b = 0, l = 0, unit = "cm")),
    axis.text.y = element_text(size = base_size, color = "black", 
                               lineheight = 0.9, 
                               margin = margin(t = 0, r = 0.05, b = 0, l = 0, unit = "cm")),
    axis.title.x = element_text(size = base_size, 
                                color = "black", 
                                margin = margin(t = 0.1, r = 0, b = 0.05, l = 0, unit = "cm")),
    # t (top), r (right), b (bottom), l (left)
    axis.title.y = element_text(size = base_size, 
                                color = "black", angle = 90, 
                                margin = margin(t = 0, r = 0.1, b = 0, l = 0.05, unit = "cm")),
    legend.background = element_rect(color = NA, 
                                     fill = NA),
    legend.key.size = unit(0.5, "lines"),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.margin = margin(t = 0, r = 0.1, b = 0, l = 0, unit = "cm"),
    legend.text = element_text(size = base_size, 
                               color = "black"),
    legend.title = element_text(size = base_size, 
                                face = "bold", 
                                hjust = 0, 
                                color = "black"),
    legend.text.align = NULL,
    legend.title.align = NULL,
    legend.direction = "vertical",
    legend.box = NULL,
    panel.background = element_rect(fill = "white", 
                                    color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = base_size * 1 ,
                              face = "bold", 
                              color = "black",
                              hjust = 0.5,
                              margin = margin(t = 0.05, b = 0, unit = "cm")),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )
}

