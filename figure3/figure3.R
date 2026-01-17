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


###figure3b-d----
#figure3b Central slice of a simulated tumor under neutral evolution, with labels indicating sample
#figure3c Maximum parsimony tree constructed from samples in figure3b
#figure3d two-dimensional distribution of spatial heterogeneity for the tumor in figure3b
tree_color <- c("#F8766D","#00BFC4","#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
work_dir <- "../data/csv_example"
state <- "Neutral"
filename <- "hu_scoef0_death0.3_adv1e-06_pushprob1_mutrate0.6_sample16vaf_cutoff0_seed54_deme_vaf.txt"
#parse simulation parameter
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
title <- paste0("fig3_hu_","scoef", s_coef, "_death", death_rate, "_adv", adv_rate, "_pushprob", push_prop, "_mutrate", mutation_rate, "_sample", sample_diameter, "vaf_cutoff", vaf_cutoff ,"_seed",  seed )
Sim_data = read.table(file = paste0(work_dir, "/", title,"_deme_vaf.txt"), header = T, stringsAsFactors = F)
boundary_strength = read_csv(file = paste0(work_dir, "/", title,"_result_boundary.csv"))
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
rpar <- readRDS(paste0(work_dir, "/", title, "_tree.rds"))
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
ggsave("figure3b.pdf", unit="cm", width=5, height=5,dpi = 300)

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
  ylim(0, max(all_point_ITH$y)+0)+
  labs(x = NULL, y = NULL, fill = "ITH") +
  theme_classic()

title_txt <- sprintf(
  "Genetic heterogeneity<br/><span style='font-size:7pt'>(%s; Silhouette = %.2f; ENN<sub>MN</sub> = 0)</span>",
  state, boundary_strength$silhouette_mean)
ITH_plot <-  p1 +
  ggnewscale::new_scale_fill()+
  geom_point(data = loc_genetic_cluster,aes(x=X,y=Y,fill= cluster),size=4,shape=21,stroke=0.4)+
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=rownames(loc_genetic_cluster)),size=5/.pt,vjust=-1.5)+
  scale_fill_manual(values=tree_color)+labs(fill="Class")+
  labs(title = title_txt) +
  guides(fill = guide_legend(order=1))+
  xlim(0,195)+
  coord_fixed()+
  my_theme()+
  theme(
    plot.title = element_markdown(lineheight = 1.3, hjust = 0.5,margin = margin(0,0,0,0)),
    legend.box.margin = margin(t = -0.2, unit = "cm"))
ggsave("figure3c.pdf", unit="cm", width=6, height=5,dpi = 300)

#plot tumor
full_location = read.table(file = paste0(work_dir, "/fig3_snapshot_dim3_pushrandom_" , title, "_prop8_8.txt"), header = T)
clone_colors = c(
  "1" = "#377EB7",
  "2" = "#E3211C")
space_plot <-  ggplot()+
  geom_tile(data=full_location, aes(x=x, y=y, fill=factor(category)),alpha=0.8, width = 1, height = 1) +
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=rownames(loc_genetic_cluster)),size=5/.pt)+
  scale_fill_manual(
    values = clone_colors,
    labels = c("1" = "Neutral", "2" = "Selected") 
  ) +
  labs(fill = "Deme", title=paste0("Tumor (", "neutral", ")"), x=NULL, y=NULL) +
  coord_fixed()+
  my_theme()
ggsave("figure3d.pdf", unit="cm", width=6, height=5,dpi = 300)

###figure3e-g----
#figure3e Central slice of a simulated tumor under selected evolution, with labels indicating sample
#figure3f Maximum parsimony tree constructed from samples in figure3e
#figure3g two-dimensional distribution of spatial heterogeneity for the tumor in figure3e
#code comment see figure3b-3d
tree_color <- c("#00BFC4","#F8766D", "#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
work_dir <- "../data/csv_example"
state <- "Selected"
filename <- "hu_scoef0.6_death0.3_adv1e-06_pushprob1_mutrate0.6_sample16vaf_cutoff0_seed31_deme_vaf.txt"
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

title <- paste0("fig3_hu_","scoef", s_coef, "_death", death_rate, "_adv", adv_rate, "_pushprob", push_prop, "_mutrate", mutation_rate, "_sample", sample_diameter, "vaf_cutoff", vaf_cutoff ,"_seed",  seed )
dir <- sprintf("scoef%.2f_death%.2f_adv%.2f_pushprob%.2f_mutrate%.2f_sample%.0f_vaf%.2f", s_coef, death_rate, adv_rate, push_prop, mutation_rate, sample_diameter, vaf_cutoff)
Sim_data = read.table(file = paste0(work_dir, "/", title,"_deme_vaf.txt"), header = T, stringsAsFactors = F)
boundary_strength = read_csv(file = paste0(work_dir, "/", title,"_result_boundary.csv"))
location = read.table(file = paste0(work_dir, "/", title,"_deme_location.txt"), header = T)
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
genetic_cluster <- get_CH_index(rpar)
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

ggsave("figure3e.pdf", unit="cm", width=5, height=5,dpi = 300)

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
  xlim(0, max(all_point_ITH$x)+10) +
  ylim(0, max(all_point_ITH$y)+0)+
  labs(x = NULL, y = NULL, fill = "ITH") +
  theme_classic()

title_txt <- sprintf(
  "Genetic heterogeneity<br/><span style='font-size:7pt'>(%s; Silhouette = %.2f; ENN<sub>MN</sub> = 0)</span>",
  state, boundary_strength$silhouette_mean)
ITH_plot <-  p1 +
  ggnewscale::new_scale_fill()+
  geom_point(data = loc_genetic_cluster,aes(x=X,y=Y,fill= cluster),size=4,shape=21,stroke=0.4)+
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=rownames(loc_genetic_cluster)),size=5/.pt,vjust=-1.5)+
  scale_fill_manual(values=tree_color)+labs(fill="Class")+
  xlim(0,160)+
  ylim(30,190)+
  labs(title = title_txt) +
  guides(fill = guide_legend(order=1))+
  coord_fixed()+
  my_theme()+
  theme(
    plot.title = element_markdown(lineheight = 1.3, hjust = 0.5,margin = margin(0,0,0,0)),
    legend.box.margin = margin(t = -0.2, unit = "cm"))
ggsave("figure3f.pdf", unit="cm", width=6, height=5,dpi = 300)

full_location = read.table(file = paste0(work_dir, "/fig3_snapshot_dim3_pushrandom_" , title, "_prop8_8.txt"), header = T)
clone_colors = c(
  "1" = "#377EB7",
  "2" = "#E3211C")
space_plot <-  ggplot()+
  geom_tile(data=full_location, aes(x=x, y=y, fill=factor(category)),alpha=0.8, width = 1, height = 1) +
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=rownames(loc_genetic_cluster)),size=5/.pt)+
  scale_fill_manual(
    values = clone_colors,
    labels = c("1" = "Neutral", "2" = "Selected")  # 重命名图例标签
  ) +
  labs(fill = "Deme", title=paste0("Tumor (", "selected", ")"), x=NULL, y=NULL) +
  coord_fixed()+
  my_theme()
ggsave("figure3g.pdf", unit="cm", width=6, height=5,dpi = 300)


###figure3h----
#Boxplot of ENNmn values across different selection intensity with varying selective coefficients
df <- read_csv("../data/fig3_hu_result_hu_time_result.csv")
df_long <- df  %>%
  pivot_longer(cols = c(value_enn, value_lsi),
               names_to = "ls_type",
               values_to = "value") %>%
  mutate(ls_type = ifelse(ls_type == "value_enn", "enn", "lsi")) %>%
  filter(ls_type == "enn")
small_df <- df_long %>% filter(death_rate==0.3 & adv_rate %in% c(0.000001,0) & push_prop==1 & mut_rate==0.6 & sample_diameter==16 & vaf_cutoff==0) %>%
  mutate(s_coef = as.character(s_coef))
comparisons <- lapply(
  unique(small_df$s_coef[small_df$s_coef != "0"]),
  function(x) c("0", x)
)
ggplot(small_df, aes(x = s_coef, y = value)) +
  geom_boxplot(aes(fill = s_coef),outliers = F) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6, color="#E64B35FF") +
  labs(x = "Selection coefficient", y = expression(ENN[MN]), fill = "s_coef", title=NULL) +
  my_theme()+theme(
    legend.position = "none")
ggsave("figure3h.pdf",unit="cm", width=4, height=5,dpi = 300)

###figure3i----
#Boxplot of Silhouette values comparing neutrally evolving versus selected tumors.
df <- read_csv("../data/fig3_hu_time_boundary_silhouette_freq.csv")
#selection
df_select <- df %>% filter(death_rate==0.3 & adv_rate %in% c(0.000001,0) & 
                             push_prop==1 & mut_rate==0.6 & sample_diameter==16 & vaf_cutoff==0 &
                             s_coef == 0.6) %>% dplyr::select(silhouette_mean, seed, cat1, cat2,whole_adv_freq) %>% 
  group_by(seed) %>% summarise(silhouette_mean=max(silhouette_mean), adv_freq= mean(whole_adv_freq))
df_select$type <- "selection"

#neutral 
df2 <- read_csv("../data/fig3_hu_time_boundary_silhouette.csv") 
df_neutral <- df2 %>% filter(death_rate==0.3 & adv_rate %in% c(0.000001,0) & 
                               push_prop==1 & mut_rate==0.6 & sample_diameter==16 & vaf_cutoff==0 &
                               s_coef == 0) %>% dplyr::select(silhouette_mean, seed, cat1, cat2) %>% 
  group_by(seed) %>% summarise(silhouette_mean=max(silhouette_mean)) 
df_neutral$type <- "neutral"
df_neutral$adv_freq <- 0

df_all <- rbind(df_neutral, df_select)

#Distinguish different selection situations according to the frequency of the selected cells
df_plot <- df_all %>% 
  mutate(
    adv_freq = 2 * adv_freq,
    adv_group = case_when(
      type == "selection" & adv_freq >= 0   & adv_freq < 0.1 ~ "0-0.1",
      type == "selection" & adv_freq >= 0.1 & adv_freq < 0.9 ~ "0.1-0.9",
      type == "selection" & adv_freq >= 0.9 ~ "0.9-1",
      type == "neutral" ~ "neutral",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(adv_group))
df_plot$adv_group <- factor(df_plot$adv_group, 
                            levels = c("neutral", "0-0.1", "0.1-0.9", "0.9-1"))

p <- ggplot(df_plot, aes(x = adv_group, y = silhouette_mean)) +
  geom_boxplot(
    outlier.shape = NA
  ) +
  geom_point(
    position = position_jitter(width = 0.15),
    size = 0.5,
    alpha = 0.6,
    color = "#E64B35FF"
  ) +
  stat_compare_means(
    comparisons = list(
      c("neutral", "0-0.1"),
      c("neutral", "0.1-0.9"),
      c("neutral", "0.9-1")
    ),
    method = "wilcox.test",
    label = "p.signif",
    size = 6/.pt
  ) +
  scale_x_discrete(labels = c(
    "neutral"  = "Neutral",
    "0-0.1"    = "Selection\n(0-0.1)",
    "0.1-0.9"  = "Selection\n(0.1-0.9)",
    "0.9-1"    = "Selection\n(0.9-1)"
  )) +
  labs(x = "", y = "Silhouette value") +
  my_theme() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(3, 0, -8, 0),
    panel.spacing = unit(0, "cm")
  )
ggsave("figure3i.pdf", units="cm", width=5, height=5)

###figure3j----
#Boxplot of ENNmn values across different pushing power
df <- read_csv("../data/fig3_hu_result_hu_time_result.csv")
small_df <- df_long %>% filter(s_coef==0 &  death_rate==0.3 & mut_rate==0.6 & sample_diameter==16 & vaf_cutoff==0 & push_prop %in% c(0, 0.5, 1)) %>%
  mutate(push_prop = as.character(push_prop))
comparisons <- lapply(
  unique(df_long$push_prop[df_long$push_prop != "0"]),
  function(x) c("0", x)
)
ggplot(small_df, aes(x = push_prop, y = value)) +
  geom_boxplot(aes(fill = push_prop), outliers = F) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6, color="#1f77b4") +
  labs(x = "Push power", y = expression(ENN[MN]), fill = "push_prop", title=NULL) +
  my_theme()+theme(  
    legend.position = "none")
ggsave("figure3j.pdf",unit="cm", width=3.5, height=5,dpi = 300)


###figure3k----
#Boxplot of Silhouette values comparing tumors with different pushing power
df <- read_csv("../data/fig3_hu_time_boundary_silhouette_push.csv")
df_all <- df %>%
  filter(
    death_rate==0.3,
    adv_rate %in% c(0, 0.000001),
    push_prop %in% c(0, 0.5, 1),
    mut_rate==0.6,
    sample_diameter==16,
    vaf_cutoff==0,
    s_coef == 0
  ) %>%
  dplyr::select(silhouette_mean, seed, cat1, cat2, push_prop)
df_all$push_prop <- factor(df_all$push_prop)
x_levels <- sort(unique(df_all$push_prop))
target <- "0"
others <- setdiff(x_levels, target)
ordered_others <- others[order(abs(as.numeric(as.character(others)) -
                                     as.numeric(target)))]

comparisons <- lapply(ordered_others, function(x) c(target, x))
ymax <- max(df_all$silhouette_mean, na.rm = TRUE)
ymin <- min(df_all$silhouette_mean, na.rm = TRUE)
step <- 0.1 * (ymax - ymin)
label_y <- ymax + step * seq_along(comparisons)

plot_all <- ggplot(df_all, aes(x = push_prop, y = silhouette_mean)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(color = "#1f77b4", width = 0.2, size = 0.5, alpha = 0.6) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    label = "p.signif",
    size = 5/.pt,
    label.y = label_y
  ) +
  labs(x = "Push power", y = "Silhouette value") +
  my_theme() +
  theme(
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(ymin, max(label_y) + step*1), clip = "off")
ggsave("figure3k.pdf", units="cm", width=3.5, height=5)