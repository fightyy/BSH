suppressPackageStartupMessages({ 
  library(tidyverse)
  library(this.path)
  library(ggtree)
  library(ggplot2)
  library(ks)
  library(spatstat)
  library(GET)
  library(phangorn)
  library(cowplot)
  library(landscapemetrics)
  library(terra)
  library(raster)
  library(patchwork)
  library(mobster)
  library(ggtext)
})

setwd(this.dir())
source("../bin/theme_setup_simple.R")
source("../bin/function.R")

###figure4a-4c----
#figure4a:Simulated (2D) tumor with migration rate (p) set to 0.1
#figure4b:The phylogenetic tree of the sampled sectors from the simulated tumor in figure4a
#figure4c:The spatial distribution of the genetic classes for the tumor in figure4a
work_dir <- "../data/csv_example"
state <- "Migration"
#parse simulation parameter
filename <- "hu_scoef0_death0.3_adv0_pushprob1_mutrate0.6_sample16_vaf_cutoff0_migrate0.1_seed74_deme_vaf.txt"
s_coef         <- str_extract(filename, "(?<=scoef)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
death_rate     <- str_extract(filename, "(?<=death)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
adv_rate       <- str_extract(filename, "(?<=adv)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
push_prop      <- str_extract(filename, "(?<=pushprob)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
mutation_rate  <- str_extract(filename, "(?<=mutrate)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
sample_diameter<- str_extract(filename, "(?<=sample)\\d+(?=[_\\.]|$)") %>% as.numeric()
vaf_cutoff     <- str_extract(filename, "(?<=cutoff)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
migrate_rate   <- str_extract(filename, "(?<=migrate)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
seed           <- str_extract(filename, "(?<=seed)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?(?=[_\\.]|$)") %>% as.numeric()
r <- (200/2)
lambda_b <- 1-death_rate
push_power <- r*push_prop
birth_rate <- (1+s_coef)*lambda_b+death_rate
tree_color <- c("#F8766D","#00BFC4","#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
title <- paste0("fig4_hu_","scoef", s_coef, "_death", death_rate, "_adv", adv_rate, "_pushprob", push_prop, "_mutrate", mutation_rate,
                "_sample", sample_diameter, "_vaf_cutoff", vaf_cutoff , "_migrate",  migrate_rate, "_seed",  seed )
dir <- sprintf("scoef%.2f_death%.2f_adv%.2f_pushprob%.2f_mutrate%.2f_sample%.0f_vaf%.2f_migrate%.2f", s_coef, death_rate, adv_rate, push_prop, mutation_rate, sample_diameter, vaf_cutoff, migrate_rate)
Sim_data = read.table(file = paste0(work_dir, "/", title,"_deme_vaf.txt"), header = T, stringsAsFactors = F)
location = read.table(file = paste0(work_dir, "/", title,"_deme_location.txt"), header = T)

#convert vaf count to vaf matrix
variant_counts = c()
depth = c()
k = 2
while (k < length(colnames(Sim_data))) {
  tmp = Sim_data[, k]*Sim_data[,k+1]
  variant_counts = cbind(variant_counts, tmp)
  depth = cbind(depth, Sim_data[, k])
  k = k + 2
}
sample_name <- location$sample
colnames(variant_counts) = sample_name
variant_counts = as.data.frame(variant_counts)
depth = as.data.frame(depth)
vaf_table <- (variant_counts/depth) %>% mutate_all(~replace(., is.nan(.), 0))
vaf_table$GL <- 0
mutated = vaf_table > vaf_cutoff
storage.mode(mutated) = "numeric"
phy_data = mutated %>% t() %>% phangorn::phyDat(type="USER", levels=c(0,1))
attr(phy_data, "id") = rownames(mutated)
rpar <- readRDS(paste0("../data/csv_example/", title, "_tree.rds"))
#Get best clustering based on CH index
genetic_cluster <- get_CH_index(rpar)
#plot tree
group <- genetic_cluster$Genetic_group
annotation <- data.frame(label=names(group),group=as.character(group),sample=names(group))
rpar1 <- drop.tip(rpar, "GL")
trs <- full_join(rpar1,annotation,by='label')
ntip <- ape::Ntip(trs)
xmax <- max(phytools::nodeHeights(rpar1))
scale_width <- 50
tree_plot <- ggtree(trs) +
  geom_tippoint(aes(color = group), shape = 16, size = 10/.pt) +
  geom_treescale(x = xmax - scale_width + 10, y = 0, width = scale_width, fontsize = 6/.pt)+
  geom_text2(aes(label = sample), size = 4/.pt, color = "white") +
  scale_color_manual(values = tree_color, name = "Genetic class") +
  geom_rootpoint() +
  labs(title = "Tree") +
  coord_cartesian(clip = "off") +
  my_theme() +
  theme(
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    axis.text.x     = element_blank(),
    axis.text.y     = element_blank(),
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank()
  )
ggsave("figure4b.pdf", unit="cm", width=4, height=5,dpi = 300)

# Map cluster labels to spatial coordinates of bulk samples
loc_genetic_cluster <- data.frame()
for(sample_ in location$sample){
  cluster <- as.character(genetic_cluster$Genetic_group[[sample_]])
  coord_vec <- location %>% filter(sample == sample_) %>% dplyr::select(x,y,z)
  x <- coord_vec[["x"]]
  y <- coord_vec[["y"]]
  z <- coord_vec[["z"]]
  small_df <- data.frame(X=x, Y=y, Z=z, cluster=cluster, sample=sample_)
  rownames(small_df) <- sample_
  loc_genetic_cluster <- rbind(loc_genetic_cluster, small_df)
}
cluster_num <- length(unique(loc_genetic_cluster$cluster))
select_info <- loc_genetic_cluster[c("X","Y")]
ITH_matrix <- get_ITH_matrix(phy_data,select_info,r)
# Random subsampling of ITH grid points weighted by ITH value
set.seed(123)
ITH <- ITH_matrix$ITH
temp_order <- c()
for (n in 1:length(ITH)){
  if (sample(c(0,1),1,prob=c(1-ITH[n],ITH[n]))==1){
    temp_order=c(temp_order,n)
  }
}
all_point_ITH <- ITH_matrix[temp_order,]

#plot ITH map
p1 <- ggplot(data=all_point_ITH,aes(x=x,y=y))+
  stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
  labs(fill="ITH")+
  scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                       breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
  xlim(0, max(all_point_ITH$x)+10) +
  ylim(0, max(all_point_ITH$y)+10)+
  labs(x = NULL, y = NULL, fill = "ITH") +
  theme_classic()
title_txt <- sprintf(
  "Genetic heterogeneity<br/><span style='font-size:7pt'>(%s; Silhouette = %.2f; ENN<sub>MN</sub> = 0.21)</span>",
  state, round(0.119145619, 2)
)
ITH_plot <-  p1 +
  ggnewscale::new_scale_fill()+
  geom_point(data = loc_genetic_cluster,aes(x=X,y=Y,fill= cluster),size=4,shape=21,stroke=0.4)+
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=rownames(loc_genetic_cluster)),size=5/.pt,vjust=-1.5)+
  scale_fill_manual(values=tree_color)+labs(fill="Class")+
  labs(title = title_txt) +
  guides(fill = guide_legend(order=1))+
  coord_fixed()+
  my_theme()+
  theme(
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    #plot.title = element_text(size = 7, hjust = 0.5,lineheight =0),
    plot.title = ggtext::element_markdown(lineheight = 1.3, hjust = 0.5,margin = margin(0,0,0,0)),
    legend.box.margin = margin(t = -0.2, unit = "cm"))
ggsave("figure4c.pdf", unit="cm", width=6, height=5,dpi = 300)

#plot tumor
full_location = read.table(file = paste0(work_dir, "/fig4_snapshot_dim2_pushrandom_" , title, "_prop8_8.txt"), header = T)
tumor_diameter <- 16L
spacing        <- 16L
step_size      <- tumor_diameter + spacing
#generate all sample point
bounds_by_x <- full_location %>%
  mutate(x = as.integer(x), y = as.integer(y)) %>%
  group_by(x) %>%
  summarise(ymin = min(y), ymax = max(y), .groups = "drop")
min_x <- min(bounds_by_x$x); max_x <- max(bounds_by_x$x)
min_y <- min(bounds_by_x$ymin); max_y <- max(bounds_by_x$ymax)
x_coords <- seq(from = min_x + tumor_diameter, to = max_x, by = step_size)
y_coords <- seq(from = min_y + tumor_diameter, to = max_y, by = step_size)
grid_pts <- expand.grid(x = as.integer(x_coords), y = as.integer(y_coords))
#Only retain the sampling points that fall within the boundary
df_sample_all <- grid_pts %>%
  inner_join(bounds_by_x, by = "x") %>%
  filter(y >= ymin, y <= ymax) %>%
  dplyr::select(x, y) %>%
  distinct()
clone_colors = c(
  "1" = "#377EB7",
  "2" = "#E3211C")
space_plot <- ggplot() +
  geom_tile(
    data = full_location,
    aes(x = x, y = y, fill = factor(category)),
    alpha = 0.8, width = 1, height = 1
  ) +
  geom_point(
    data = df_sample_all,
    aes(x = x, y = y, color = "low"),
    size = 1.8, stroke = 0
  ) +
  geom_point(
    data = loc_genetic_cluster,
    aes(x = X, y = Y, color = "high"),
    size = 1.8, stroke = 0
  ) +
  geom_text(
    data = loc_genetic_cluster,
    aes(x = X, y = Y, label = rownames(loc_genetic_cluster)),
    size = 6/.pt, vjust = -1, color = "#E64B35FF"
  ) +
  scale_fill_manual(
    values = clone_colors,
    labels = c("1" = "Neutral", "2" = "Selection"),
    name   = "Deme"
  ) +
  scale_color_manual(
    name   = "Deme density",
    values = c(high = "#E64B35FF", low = "grey50"),
    breaks = c("high","low"),
    labels = c("high" = "High", "low" = "Low")
  ) +
  guides(
    fill  = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  labs(title = "Tumor\n(Migration rate = 0.1)", x = NULL, y = NULL) +
  coord_fixed()+
  my_theme()
ggsave("figure4a.pdf", unit="cm", width=6, height=5,dpi = 300)


###figure4d----
#The ENNmn values for simulated tumors with different migration rates
df <- read_csv("../data/fig4_result_enn_treestat_migrate.csv")
df_long <- df  %>%
  pivot_longer(cols = c(value_enn, value_lsi),
               names_to = "ls_type",
               values_to = "value") %>%
  mutate(ls_type = ifelse(ls_type == "value_enn", "enn", "lsi")) %>%
  filter(ls_type == "enn")
small_df <- df_long %>%
  filter(
    death_rate == 0.3,
    adv_rate %in% c(0.000001, 0),
    push_prop == 1,
    mut_rate == 0.6,
    sample_diameter == 16,
    vaf_cutoff == 0,
    s_coef == 0
  ) %>%
  mutate(migrate_rate_num = as.numeric(migrate_rate))
x_levels <- sort(unique(small_df$migrate_rate_num))
small_df <- small_df %>%
  mutate(migrate_rate = factor(migrate_rate_num, levels = x_levels, labels = as.character(x_levels)))
target <- 0.5
others <- setdiff(x_levels, target)
ordered_others <- others[order(abs(others - target))]
comparisons <- lapply(ordered_others, function(x) c(as.character(target), as.character(x)))
y_max   <- max(small_df$value, na.rm = TRUE)
y_min   <- min(small_df$value, na.rm = TRUE)
step    <- 0.15 * (y_max - y_min)
label_y <- y_max + step * seq_along(comparisons)
#plot
p <- ggplot(small_df, aes(x = migrate_rate, y = value)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6, color = "#ff7f0e") +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    label = "p.signif",
    size = 7/.pt,
    label.y = label_y
  ) +
  labs(
    x = "Migration rate",
    y = expression(ENN[MN]),
    title = "Deme model (2D simulation)"
  ) +
  my_theme() +
  theme(
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y  = element_text(size = 7),
    axis.title.y = element_text(size = 9),
    strip.text   = element_text(size = 7, margin = margin(t = 1, b = 1)),
    legend.position = "none",
    plot.title = element_text(hjust=1.2),
    plot.margin = margin(0.1, 0, 0, 0, "cm")
  ) +
  coord_cartesian(ylim = c(y_min, max(label_y) + step))
ggsave("figure4d.pdf", units="cm", width=3.5, height=5)


##figure4e----
#The Silhouette values for simulated tumors with different migration rates
df <- read_csv("../data/fig4_hu_time_boundary_silhouette_migration.csv")
small_df <- df %>%
  filter(
    death_rate == 0.3,
    adv_rate %in% c(0.000001, 0),
    push_prop == 1,
    mut_rate == 0.6,
    sample_diameter == 16,
    vaf_cutoff == 0,
    s_coef == 0,
    migrate_rate %in% c("0", "0.5" )
  ) %>%
  mutate(migrate_rate_num = as.numeric(migrate_rate))
x_levels <- sort(unique(small_df$migrate_rate_num))
small_df <- small_df %>%
  mutate(migrate_rate = factor(migrate_rate_num, levels = x_levels, labels = as.character(x_levels)))
x_levels_f <- levels(small_df$migrate_rate)
comparisons <- combn(x_levels_f, 2, simplify = FALSE)
p <- ggplot(small_df, aes(x = migrate_rate, y = silhouette_mean)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6, color = "#ff7f0e") +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    label = "p.signif",
    size = 7/.pt
  ) +
  labs(
    x = "Migration rate",
    y = "Silhouette value",
    title = "Deme model (2D simulation)"
  ) +
  ylim(0, 0.4)+
  my_theme() +
  theme(
    panel.spacing.x = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.y  = element_text(size = 7),
    axis.title.y = element_text(size = 9),
    strip.text   = element_text(size = 7, margin = margin(t = 1, b = 1)),
    plot.title = element_text(hjust=1.2),
    legend.position = "none",
    plot.margin = margin(0.1, 0, 0, 0, "cm")
  )
ggsave("figure4e.pdf", units="cm", width=3.5, height=5)

###figure4f----
#GO enrichment of the differentially expressed genes between HCC tumors with high ENNmn values (weak blockness) and low ENNmn values (strong blockness)
df_enrich <- read_csv("../data/fig4_go_enrich.csv") %>% mutate(score = abs(score))
df_enrich_sorted <- df_enrich[order(df_enrich$score, decreasing = FALSE), ]
top5_up <- df_enrich[df_enrich$direction == "Up", ] |>
  dplyr::arrange(pvalue) |>
  head(3)
top5_down <- df_enrich[df_enrich$direction == "Down", ] |>
  dplyr::arrange(pvalue) |>
  head(3)
df_enrich_plot <-  rbind(top5_down, top5_up) %>%
  mutate(Description_wrap = str_wrap(Description, width = 39)) %>%
  group_by(direction) %>%
  mutate(
    Description_wrap = fct_reorder(
      Description_wrap,
      score
    )
  ) %>%
  ungroup()
title_txt <- sprintf(
  "GO enrichment for HCC (high ENN<sub>MN</sub> VS low ENN<sub>MN</sub>)")
ggplot(df_enrich_plot,
       aes(x = score, y = Description_wrap, fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c("Down" = "#1f77b4",
               "Up"   = "#d62728")
  ) +
  scale_x_continuous(labels = scales::label_number())+
  facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Enrichment significance (−log10 (q))",
    y = NULL,
    fill = NULL,
    title = title_txt
  ) +
  my_theme()+
  theme(
    strip.text.y = element_text(face = "bold", size =6 ),
    plot.title = element_markdown(lineheight = 1.3, hjust = 0.8 ,margin = margin(0,0,0,0)),
    axis.text.y  = element_text(size = 4),
    legend.position = "none"
  )
ggsave("figure4f.pdf", unit="cm", width=6.5, height=5,dpi = 300)

###figure4g
#Regression between GSVA score of the GO:0098742 (“cell–cell adhesion via plasma-membrane adhesion molecules”) vs the ENNmn values
df_score_enn <- read_csv("../data/fig4_gsva_enn.csv")
fit <- lm(GO_0098742_score ~ value_enn, data = df_score_enn)
pval <- summary(fit)$coefficients[2, 4]

#P-value in scientific notation
decompose_p <- function(p) {
  if (is.na(p) || p <= 0) {
    p <- .Machine$double.xmin
  }
  k <- floor(log10(p))
  mant <- p / (10^k)
  list(mant = mant, expo = k)
}
dp <- decompose_p(pval)
mant_str <- formatC(dp$mant, format = "f", digits = 2)
text_label_expr <- bquote(
  "p" == .(mant_str) %*% 10^.(dp$expo) )
text_label_str <- as.character(as.expression(text_label_expr))

#plot
ggplot(df_score_enn, aes(x = value_enn, y = GO_0098742_score)) +
  geom_point(color="#1f77b4",size = 1, alpha = 0.5) +
  geom_smooth(method = "lm", se = T,linewidth = 0.8, color = "red", linetype = 2) +
  annotate("text", x = mean(range(df_score_enn$value_enn)), y = max(df_score_enn$GO_0098742_score),
           label = text_label_str, parse = TRUE,
           size = 7/.pt, fontface = "plain", hjust = 0.5) +
  labs(
    x = expression(ENN[MN]),
    y = "GSVA score (GO:0098742)",
    title = "HCC"
  )+
  my_theme()+
  theme(
    legend.position = "none")
ggsave("figure4g.pdf", width=5, height=5, units="cm", dpi = 300)



###figure4h-4j----
#figure4h:A neutrally evolving tumor simulated using the boundary growth model
#figure4i:The phylogenetic tree of the sampled sectors from the simulation in figure4h
#figure4j:The spatial distribution of the genetic classes for the tumor in panel figure4h
#code comment see figure4a-4c
tree_color <- c("#F8766D","#00BFC4","#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
work_dir <- "../data/csv_example"
filename <- "hu_scoef0_death0.3_adv0_pushprob0_mutrate0.6_sample16vaf_cutoff0_seed5_deme_vaf.txt"
state <- "Boundary growth"
s_coef         <- str_extract(filename, "(?<=scoef)\\d+(\\.\\d+)?") %>% as.numeric()
death_rate     <- str_extract(filename, "(?<=death)\\d+(\\.\\d+)?") %>% as.numeric()
adv_rate <- str_extract(filename, "(?<=adv)[+-]?\\d+(\\.\\d+)?([eE][+-]?\\d+)?") %>% as.numeric()
push_prop      <- str_extract(filename, "(?<=pushprob)\\d+(\\.\\d+)?") %>% as.numeric()
mutation_rate  <- str_extract(filename, "(?<=mutrate)\\d+(\\.\\d+)?") %>% as.numeric()
sample_diameter<- str_extract(filename, "(?<=sample)\\d+(?=vaf)") %>% as.numeric()
vaf_cutoff     <- str_extract(filename, "(?<=cutoff)\\d+(\\.\\d+)?") %>% as.numeric()
seed     <- str_extract(filename, "(?<=seed)\\d+(\\.\\d+)?") %>% as.numeric()
r <- (200/2)
lambda_b <- 1-death_rate
push_power <- r*push_prop
birth_rate <- (1+s_coef)*lambda_b+death_rate
title <- paste0("fig4_hu_","scoef", s_coef, "_death", death_rate, "_adv", adv_rate, "_pushprob", push_prop, "_mutrate", mutation_rate, "_sample", sample_diameter, "vaf_cutoff", vaf_cutoff ,"_seed",  seed )
dir <- sprintf("scoef%.2f_death%.2f_adv%.2f_pushprob%.2f_mutrate%.2f_sample%.0f_vaf%.2f", s_coef, death_rate, adv_rate, push_prop, mutation_rate, sample_diameter, vaf_cutoff)
Sim_data = read.table(file = paste0(work_dir, "/", title,"_deme_vaf.txt"), header = T, stringsAsFactors = F)
boundary_strength = read_csv(file = paste0(work_dir, "/", title,"_result_boundary.csv"))
location = read.table(file = paste0(work_dir, "/", title,"_deme_location.txt"), header = T)
all_location = read.table(file = paste0(work_dir, "/", "fig4_snapshot_dim3_pushrandom_", title,"_prop8_8.txt"), header = T)
x_rng <- range(all_location$x)
y_rng <- range(all_location$y)
max_distance <- sqrt((x_rng[2] - x_rng[1])^2 + (y_rng[2] - y_rng[1])^2)
boundary <- all_location %>%
  group_by(y) %>%
  summarise(
    x_min = min(x),
    x_max = max(x),
    .groups = "drop"
  )
boundary_pts <- bind_rows(
  boundary %>% transmute(x = x_min, y),
  boundary %>% transmute(x = x_max, y)
)
#Calculate the minimum distance from each sample to all boundary points
sample_pts <- location %>% dplyr::select(sample, x, y)

min_distances <- sample_pts %>%
  rowwise() %>%
  mutate(
    min_distance = min(sqrt((x - boundary_pts$x)^2 + (y - boundary_pts$y)^2))
  ) %>%
  ungroup()

annotated <- min_distances %>%
  mutate(
    location = ifelse(min_distance > 0.25 * max_distance, "center", "edge"),
    max_distance = max_distance)

#plot tree
#convert vaf count to vaf matrix
variant_counts = c()
depth = c()
k = 2
while (k < length(colnames(Sim_data))) {
  tmp = Sim_data[, k]*Sim_data[,k+1]
  variant_counts = cbind(variant_counts, tmp)
  depth = cbind(depth, Sim_data[, k])
  k = k + 2
}
sample_name <- location$sample
colnames(variant_counts) = sample_name
variant_counts = as.data.frame(variant_counts)
depth = as.data.frame(depth)
vaf_table <- (variant_counts/depth) %>% mutate_all(~replace(., is.nan(.), 0))
vaf_table$GL <- 0
mutated = vaf_table > 0.05
storage.mode(mutated) = "numeric"
phy_data = mutated %>% t() %>% phangorn::phyDat(type="USER", levels=c(0,1))
attr(phy_data, "id") = rownames(mutated)
rpar <- readRDS(paste0("../data/csv_example/", title, "_tree.rds"))
genetic_cluster <- get_CH_index(rpar)
group <- genetic_cluster$Genetic_group
annotation <- data.frame(label=names(group),group=as.character(group),sample=names(group))
rpar1 <- drop.tip(rpar, "GL")
annotation2 <- annotation %>%
  left_join(annotated %>% dplyr::select(sample, location), by = "sample")
xmax <- max(phytools::nodeHeights(rpar1))
scale_width <- 150
tree_plot <- ggtree(rpar1) %<+% annotation2 +
  geom_tippoint(aes(color = location), shape = 16, size = 10/.pt) +
  geom_treescale(x = xmax - scale_width + 10, y = 0, width = scale_width, fontsize = 6/.pt)+
  geom_text2(aes(label = sample), size = 4/.pt, color = "white") +
  scale_color_manual(values = c(center = "#E64B35FF", edge = "#1f77b4"),
                     labels=c("center" = "Center", "edge" = "Edge"),
                     name = "Location") +
  geom_rootedge(rootedge = rpar$edge.length[length(rpar$edge.length)]) +
  geom_rootpoint() +
  coord_cartesian(clip = "off") +
  labs(title = "Tree") +
  my_theme() +
  theme(
    axis.line       = element_blank(),
    axis.ticks      = element_blank(),
    axis.text.x     = element_blank(),
    axis.text.y     = element_blank(),
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    legend.position = "right",
    legend.box.margin = margin(t = -0.5, unit = "cm")
  )
ggsave("figure4i.pdf", unit="cm", width=4, height=5,dpi = 300)

#plot ITH
loc_genetic_cluster <- data.frame()
for(sample_ in location$sample){
  cluster <- as.character(genetic_cluster$Genetic_group[[sample_]])
  coord_vec <- location %>% filter(sample == sample_) %>% dplyr::select(x,y,z)
  x <- coord_vec[["x"]]
  y <- coord_vec[["y"]]
  z <- coord_vec[["z"]]
  small_df <- data.frame(X=x, Y=y, Z=z, cluster=cluster, sample=sample_)
  rownames(small_df) <- sample_
  loc_genetic_cluster <- rbind(loc_genetic_cluster, small_df)
}
cluster_num <- length(unique(loc_genetic_cluster$cluster))
select_info <- loc_genetic_cluster[c("X","Y")]
ITH_matrix <- get_ITH_matrix(phy_data,select_info,r)
set.seed(123)
ITH <- ITH_matrix$ITH
temp_order <- c()
for (n in 1:length(ITH)){
  if (sample(c(0,1),1,prob=c(1-ITH[n],ITH[n]))==1){
    temp_order=c(temp_order,n)
  }
}
all_point_ITH <- ITH_matrix[temp_order,]
p1 <- ggplot(data=all_point_ITH,aes(x=x,y=y))+
  stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
  labs(fill="ITH")+
  scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                       breaks=c(0,0.25,0.5,0.75,1),colours =  colormap)+
  xlim(0, 2*r) +
  ylim(0, 2*r)+
  labs(x = NULL, y = NULL, fill = "ITH") +
  theme_classic()
title_txt <- sprintf(
  "Genetic heterogeneity<br/><span style='font-size:7pt'>(%s; Silhouette = %.2f; ENN<sub>MN</sub> = 0)</span>",
  state, round(max(boundary_strength$silhouette_mean), 2))
ITH_plot <-  p1 +
  ggnewscale::new_scale_fill()+
  geom_point(data = loc_genetic_cluster,aes(x=X,y=Y,fill= cluster),size=4,shape=21,stroke=0.4)+
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=rownames(loc_genetic_cluster)),size=5/.pt,vjust=-1.5)+
  scale_fill_manual(values=tree_color)+labs(fill="Class")+
  labs(title = title_txt) +
  guides(fill = guide_legend(order=1))+
  coord_fixed()+
  my_theme()+
  theme(
    plot.title = element_markdown(lineheight = 1.3, hjust = 0.4,margin = margin(0,0,0,0)),
    legend.box.margin = margin(t = -0.2, unit = "cm"),
    plot.margin = margin(t = 0.1, unit = 'cm'))
ggsave("figure4j.pdf", unit="cm", width=6, height=5,dpi = 300)

#plot tumor
full_location = read.table(file = paste0("../data/csv_example/fig4_snapshot_dim3_pushrandom_" , title, "_prop8_8.txt"), header = T)
clone_colors = c(
  "1" = "#377EB7",
  "2" = "#E3211C")
space_plot <- ggplot() +
  geom_tile(
    data = full_location,
    aes(x = x, y = y, fill = mutation_count),
    width = 1, height = 1
  ) +
  geom_text(
    data = loc_genetic_cluster,
    aes(x = X, y = Y, label = rownames(loc_genetic_cluster)),
    size = 5/.pt
  ) +
  scale_fill_gradient(
    low = "white", high = "#377EB7",
    name = "Mutation count\n(Deme level)"
  ) +
  labs(title = paste0("Tumor (boundary growth)"), x = NULL, y = NULL) +
  coord_fixed() +
  my_theme()
ggsave("figure4h.pdf", unit="cm", width=6, height=5,dpi = 300)



#figure4k----
#the Mutation ratio (edge vs center) for the simulated data as well as the real tumor.
df_boxplot_all <- read_csv("../data/fig4_result_center_edge_real_simulation.csv") %>% 
                  mutate( tumor_type = factor( tumor_type,levels = c("Simulation", "HCC", "LUAD", "CRC", "GOA") ))
lv <- levels(df_boxplot_all$tumor_type)
sim_idx <- which(lv == "Simulation")
idx_all <- seq_along(lv)
idx_others <- setdiff(idx_all, sim_idx)
ordered_others <- idx_others[order(abs(idx_others - sim_idx))]
comparisons_list <- lapply(lv[ordered_others], function(x) c("Simulation", x))
ymin <- min(df_boxplot_all$edge_center_ratio, na.rm = TRUE)
ymax <- max(df_boxplot_all$edge_center_ratio, na.rm = TRUE)
step <- 0.1 * (ymax - ymin)
label_y <- ymax + step * seq_along(comparisons_list)

p <- ggplot(df_boxplot_all, aes(x = tumor_type, y = edge_center_ratio, fill = tumor_type)) +
  geom_boxplot(outliers = FALSE) +
  geom_hline(yintercept = 1, linetype = 2, color = "red") +
  scale_fill_manual(values = c(
    "Simulation" = "#999999",
    "HCC" = "#E64B35FF",
    "LUAD" = "#1f77b4",
    "CRC" = "#2ca02c",
    "GOA" = "#ff7f0e",
    "ESCA" = "#FBB4AE",
    "BLCA" = "#B3CDE3",
    "PCa" = "#CCEBC5"
  )) +
  stat_compare_means(
    comparisons = comparisons_list,
    method = "wilcox.test",
    paired = FALSE,
    label = "p.signif",
    size = 6/.pt,
    label.y = label_y
  ) +
  labs(x = NULL, y = "Mutation ratio (edge vs center)") +
  my_theme() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1,size=5)
  ) +
  coord_cartesian(ylim = c(0.7, max(label_y) + step*1))

ggsave("figure4k.pdf",unit="cm", width=5, height=5,dpi = 300)

#figure4l----
#The ITH ratio between the edge vs central areas for the simulated tumors as well as real tumors
df <- read_csv("../data/fig4_real_data_center_edge_ITH.csv")

#real
tumor_dict <- setNames(
  c("LUAD", "HCC", "CRC", "GC", "GOA", "BLCA", "Panc-NEC"),
  c("lung_cancer", "liver_cancer", "colorectal_cancer",
    "gastric_cancer", "gastroesophageal_cancer",
    "bladder_cancer", "pancreatic_cancer")
)
edge_center_ratio_df <- df %>%
  pivot_wider(names_from = location, values_from = ITH_mean) %>%
  mutate(edge_center_ratio = edge / center) %>% na.omit() %>% 
  mutate(tumor_type = tumor_dict[tumor_type]) %>% 
  dplyr::filter(tumor_type %in% c("HCC", "LUAD", "CRC", "GOA"))
df_boxplot <- edge_center_ratio_df %>% mutate(
  tumor_type = factor(
    tumor_type,
    levels = c("HCC", "LUAD", "CRC", "GOA", "ESCA", "BLCA", "PCa", "Panc-NEC",  "GC")
  ))  %>% dplyr::select(edge_center_ratio, patient, tumor_type)

#simulation
df_sim <- read_csv("../data/fig4_simulation_center_edge_ITH.csv") %>% 
  mutate(edge_center_ratio =  ITH_edge / ITH_center) %>% na.omit() 
df_sim_boxplot <- df_sim %>% mutate(tumor_type="Simulation") %>% dplyr::rename("patient"="seed")
df_sim_boxplot$patient <- as.character(df_sim_boxplot$patient)
df_sim_boxplot <- df_sim_boxplot %>% dplyr::select(edge_center_ratio, patient, tumor_type)

#combine real and simulted data
df_boxplot_all <- rbind(df_boxplot, df_sim_boxplot)
df_boxplot_all <- df_boxplot_all %>% mutate(
  tumor_type = factor(
    tumor_type,
    levels = c("Simulation", "HCC", "LUAD", "CRC", "GOA", "ESCA", "BLCA", "PCa", "Panc-NEC",  "GC")
  )) %>% filter(!tumor_type %in% c("Panc-NEC","GC"))
all_types <- c("HCC", "LUAD", "CRC", "GOA")
comparisons_list <- lapply(setdiff(all_types, "Simulation"), function(x) c("Simulation", x))

#calculate p value 
p_dat <- ggpubr::compare_means(
  edge_center_ratio ~ tumor_type,
  data    = df_boxplot_all,
  method  = "wilcox.test",
  paired  = FALSE,
  ref.group = "Simulation"
) %>%
  dplyr::filter(group2 %in% vapply(comparisons_list, `[[`, "", 2)) %>%
  dplyr::arrange(match(group2, vapply(comparisons_list, `[[`, "", 2))) %>%
  # 手动转成星号
  dplyr::mutate(
    label = cut(p.adj,
                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                labels = c("****", "***", "**", "*", "ns"))
  )

order_g2 <- vapply(comparisons_list, `[[`, "", 2)
p_dat <- p_dat %>% dplyr::arrange(match(group2, order_g2))
fixed_y <- seq(1.35, 2.7, 0.15)   
p_dat$y.position <- fixed_y[seq_len(nrow(p_dat))]

#plot
p <- ggplot(df_boxplot_all, aes(x = tumor_type, y = edge_center_ratio, fill = tumor_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = 2, color = "red") +
  scale_fill_manual(values = c(
    "Simulation"="#999999","HCC"="#E64B35FF","LUAD"="#1f77b4","CRC"="#2ca02c",
    "GOA"="#ff7f0e","ESCA"="#FBB4AE","BLCA"="#B3CDE3","PCa"="#CCEBC5",
    "Panc-NEC"="#DECBE4","GC"="#FED9A6"
  )) +
  stat_pvalue_manual(
    data = p_dat,
    label = "label",
    xmin  = "group1",
    xmax  = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.3,
    inherit.aes = FALSE,
    size = 6 / .pt,
  ) +
  labs(x = NULL, y = "ITH ratio (edge vs center)") +
  my_theme() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size=5)) +
  coord_cartesian(ylim = c( min(df_boxplot_all$edge_center_ratio, na.rm = TRUE) - 0.05,
                            max(p_dat$y.position) + 0.1 ),
                  clip = "off") 
ggsave("figure4l.pdf",unit="cm", width=5, height=5,dpi = 300)
