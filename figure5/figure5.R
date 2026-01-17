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
  library(purrr)
  library(spdep)
  library(ggpubr)
})
#source("~/yy_application/R/bin/theme_setup.R")
setwd(this.dir())
source("../bin/theme_setup_simple.R")
source("../bin/function.R")

##figure5a-b----
#figure5a:A simulated tumor with an ongoing natural selection using the volume growth model
#figure5b:The spatial distribution of the genetic classes for the tumor in panel 5a a
tree_color <- c( "#F8766D", "#00BFC4",  "#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
work_dir <- "../data/csv_example"
#parse simulation parameter
filename <- "hu_scoef0.4_death0.3_adv1e-06_pushprob1_mutrate0.6_sample8vaf_cutoff0_seed94_deme_vaf.txt"
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
title <- paste0("fig5_hu_","scoef", s_coef, "_death", death_rate, "_adv", adv_rate, "_pushprob", push_prop, "_mutrate", mutation_rate, "_sample", sample_diameter, "vaf_cutoff", vaf_cutoff ,"_seed",  seed )
dir <- sprintf("scoef%.2f_death%.2f_adv%.2f_pushprob%.2f_mutrate%.2f_sample%.0f_vaf%.2f", s_coef, death_rate, adv_rate, push_prop, mutation_rate, sample_diameter, vaf_cutoff)
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

rpar = readRDS(paste0("../data/csv_example/", title, "_tree.rds"))

#Get best clustering based on CH index获取最佳分组
genetic_cluster <- get_CH_index(rpar)

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
  labs(x = "X", y = "Y", fill = "ITH") +
  theme_classic()
ITH_plot <-  p1 +
  ggnewscale::new_scale_fill()+
  geom_point(data = loc_genetic_cluster,aes(x=X,y=Y,fill= cluster),size=3,shape=21,stroke=0.4)+
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=rownames(loc_genetic_cluster)),size=5/.pt,vjust=-1.5)+
  scale_fill_manual(values=tree_color)+labs(fill="Class", x=NULL, y=NULL)+
  labs(title="Genetic heterogeneity (selected)" )+
  guides(fill = guide_legend(order=1))+
  my_theme()
ggsave("figure5b.pdf", unit="cm", width=6, height=5,dpi = 300)

#plot tumor
full_location = read.table(file = paste0("../data/csv_example/fig5_snapshot_dim3_pushrandom_" , title, "_prop8_8.txt"), header = T)
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
  labs(fill = "Deme", title="Tumor (selected)",x=NULL, y=NULL) +
  coord_fixed()+
  my_theme()+
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t=-0.2,unit="cm")
  )
ggsave("figure5a.pdf", unit="cm", width=5, height=5,dpi = 300)



###figure5c----
#The result from the MOBSTER run
mobster_df_single <- read_csv("../data/csv_example/fig5_hu_scoef0.4_death0.3_adv1e-06_pushprob1_mutrate0.6_sample8vaf_cutoff0_seed94_mobster_result_single.csv")
mobster_plot <- ggplot(mobster_df_single, aes(x = x, y = y)) +
  geom_point(data=mobster_df_single, aes(fill=mobster_result_vec) , size = 2, shape=21, color = "black") +
  geom_text(data =mobster_df_single,aes(x=x,y=y,label=sample),size=5/.pt,vjust=-1.2)+
  scale_fill_manual(values = c(
    "neutral" = "lightblue",
    "selected" = "red"),
    labels = c("selected" = "Selected", "neutral"="Neutral")
  ) +
  xlim(15,170)+
  ylim(30,200)+
  coord_fixed() +
  labs(
    title = "Single sample (selected)",
    fill = "MOBSTER result",
    x= NULL,
    y=NULL )+my_theme()+
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = -0.3, unit = "cm")
  )
ggsave("figure5c.pdf", unit="cm", width=6, height=5,dpi = 300)

###figure5d----
#Proportion of sectors that are tested positive for natural selection among sectors that are at the block boundary  vs those at the block-interior 
df <- read_csv("../data/fig5_hu_time_boundary_mob_selection_single_birth0.4_0.1-0.9.csv")
#find genetic class boundary according to Gabriel neighborhood graph (Samples that are neighbors to each other belong to different genetic class)
tol_value <- 1  #tolerance value (tolerate location errors)
df_boundary <- df %>%
  dplyr::filter(!is.na(x), !is.na(y)) %>%
  dplyr::mutate(sample = as.character(sample)) %>%
  dplyr::group_by(seed) %>%     
  dplyr::group_modify(~{
    dat <- .x
    n   <- nrow(dat)
    if (n < 2) return(dplyr::mutate(dat, boundary = FALSE))

    coords <- as.matrix(dat[, c("x","y")])
    rownames(coords) <- seq_len(n)

    nb <- build_gabriel_nb_tol(coords, tol = tol_value, fallback_n = 400)

    # boundary:any neighbor sample_type_vec is different from itself -> TRUE
    boundary <- vapply(seq_along(nb), function(i) {
      nei <- nb[[i]]
      if (length(nei) == 0) return(FALSE)
      any(dat$sample_type_vec[nei] != dat$sample_type_vec[i])
    }, logical(1))

    dat$boundary <- boundary
    dat
  }) %>%
  dplyr::ungroup()

#summary seleted sample proportion
df_summary <- df_boundary %>%
  group_by(seed, boundary) %>%
  summarise(selected_prop = mean(mobster_result_vec == "selected"),
            .groups = "drop") %>% mutate(boundary=ifelse( boundary ==TRUE, "block boundary", "block interior"))

p <- ggplot(df_summary, aes(x = boundary, y = selected_prop, fill=boundary)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Proportion of \nselected region detected by MOBSTER", title="Single sample (selected)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",   
                     label.y = max(df_summary$selected_prop, na.rm = TRUE) * 0.8,
                     comparisons = list(c("block boundary", "block interior")),
                     size = 7 / .pt)+
  scale_fill_manual(values = c("block boundary" = "#E69F00", "block interior" = "#56B4E9")) +
  scale_x_discrete(labels=c("block boundary" = "Block\nboundary", "block interior" = "Block\ninterior"))+
  ylim(0,0.7)+
  my_theme()+
  theme(
    legend.position = "none",
  )
ggsave("figure5d.pdf",unit="cm", width=4, height=5,dpi = 300)






#figure5e----
#Proportion of sectors that are tested positive for natural selection (y-axis) among sectors that are at the block boundary vs those at the block-interior
mobster_df <- read_csv("../data/csv_example/fig5_hu_scoef0.4_death0.3_adv1e-06_pushprob1_mutrate0.6_sample8vaf_cutoff0_seed94_mobster_result_window.csv")
mobster_df <- mobster_df %>% mutate(sample_type_vec = ifelse(sample_type_vec=="mixed", "mixed", "pure"))
mobster_plot <- ggplot(mobster_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = mobster_result_vec), color = "white", size = 0.5) +  
  geom_point(aes(shape = sample_type_vec), size = 2, color = "black") + 
  scale_fill_manual(values = c(
    "neutral" = "lightblue",
    "selected" = "red"
  ),
  labels = c("selected" = "Selected", "neutral"="Neutral")) +
  scale_shape_manual(values = c(
    "mixed" = 17, 
    "pure" = 16),
    labels = c("mixed" = "Mixed", "pure"="Pure")) +
  coord_fixed() +
  labs(
    title = "Sliding window (selected)",
    fill = "MOBSTER\nResult",
    shape = "Sample\nType",
    x=NULL,
    y=NULL
  )+my_theme()

mobster_plot
ggsave("figure5e.pdf", unit="cm", width=6.5, height=5,dpi = 300)

###figure5f----
#Proportion of windows (i.e. composite samples) that are tested positive for natural selection.
df <- read_csv("../data/fig5_hu_time_boundary_mob_selection_window_birth0.4_0.1-0.9.csv")
df2 <- df %>% group_by(seed) %>%
  summarise(mixed_select = sum(sample_type_vec == "mixed" & mobster_result_vec=="selected" ),
            mixed_neutral = sum(sample_type_vec == "mixed" & mobster_result_vec=="neutral" ),
            pure_select = sum(sample_type_vec != "mixed" & mobster_result_vec=="selected" ),
            pure_neutral = sum(sample_type_vec != "mixed" & mobster_result_vec=="neutral" )) %>%
  mutate(
    mixed_denom = mixed_select + mixed_neutral,
    pure_denom  = pure_select  + pure_neutral,
    mixed_prop  = ifelse(mixed_denom > 0, mixed_select / mixed_denom, NA_real_),
    pure_prop   = ifelse(pure_denom  > 0, pure_select  / pure_denom,  NA_real_)
  ) %>%
  dplyr::select(seed, mixed_prop, pure_prop)

plot_df <- df2 %>%
  pivot_longer(
    cols      = c(mixed_prop, pure_prop),
    names_to  = "group",    
    values_to = "proportion"
  ) %>%
  mutate(
    group = case_when(
      group == "mixed_prop" ~ "mixed",
      group == "pure_prop"  ~ "pure"
    )
  ) %>%
  filter(!is.na(proportion))

p <- ggplot(plot_df, aes(x = group, y = proportion, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape=NA) +
  scale_fill_manual(values = c("mixed" = "#E69F00", "pure" = "#56B4E9")) +
  scale_x_discrete(labels=c("mixed"="Mixed", "pure"="Pure"))+
  stat_compare_means(
    label = "p.signif",
    method     = "wilcox.test",
    paired     = FALSE,
    comparisons = list(c("mixed", "pure")),
    label.y    = max(plot_df$proportion, na.rm = TRUE) + 0.01,
    size = 7/.pt
  ) +
  labs(
    title = "Sliding window (selected)",
    x     = NULL,
    y     = "Proportion of windows\ndetected as under selection by MOBSTER"
  ) +
  ylim(0,1.1)+
  my_theme()+
  theme(
    legend.position = "none"
  )
ggsave("figure5f.pdf", unit="cm", width=4, height=5,dpi = 300)

##figure5g----
df <- read_csv("../data/fig5_mob_single_window_neutral.csv") %>% filter(mobster_test_type == "single")
#find genetic class boundary according to Gabriel neighborhood graph (Samples that are neighbors to each other belong to different genetic class)
tol_value <- 1  
df_boundary <- df %>%
  dplyr::filter(!is.na(x), !is.na(y)) %>%
  dplyr::mutate(sample = as.character(sample)) %>%
  dplyr::group_by(seed) %>%
  dplyr::group_modify(~{
    dat <- .x
    n   <- nrow(dat)
    if (n < 2) return(dplyr::mutate(dat, boundary = FALSE))
    coords <- as.matrix(dat[, c("x","y")])
    rownames(coords) <- seq_len(n)
    nb <- build_gabriel_nb_tol(coords, tol = tol_value, fallback_n = 400)
    # boundary:any neighbor sample_type_vec is different from itself -> TRUE
    boundary <- vapply(seq_along(nb), function(i) {
      nei <- nb[[i]]
      if (length(nei) == 0) return(FALSE)
      any(dat$sample_type_vec[nei] != dat$sample_type_vec[i])
    }, logical(1))

    dat$boundary <- boundary
    dat
  }) %>%
  dplyr::ungroup()

df_summary <- df_boundary %>%
  group_by(seed, boundary) %>%
  summarise(selected_prop = mean(mobster_result_vec == "selected"),
            .groups = "drop") %>% mutate(boundary=ifelse(boundary ==TRUE, "block boundary", "block interior"))

p <- ggplot(df_summary, aes(x = boundary, y = selected_prop, fill=boundary)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Proportion of \nselected region detected by MOBSTER", title="Single sample (neutral)") +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif", 
                     label.y = max(df_summary$selected_prop, na.rm = TRUE) * 1.05,
                     comparisons = list(c("block boundary", "block interior")),
                     size = 7 / .pt)+
  scale_fill_manual(values = c("block boundary" = "#E69F00", "block interior" = "#56B4E9")) +
  scale_x_discrete(
    labels = c("block boundary" = "Block\nboundary", "block interior" = "Block\ninterior")
  ) +
  ylim(0,0.4)+
  my_theme()+
  theme(
    legend.position = "none")
ggsave("figure5g.pdf",unit="cm", width=5, height=5,dpi = 300)


###figure5h----
#Proportion of sectors that are tested positive for natural selection (y-axis) when simulating neutrally evolving tumors.
#Comparison between sectors that are at the tumor periphery and those that are the tumor interior
df <- read_csv("../data/fig5_mob_single_window_neutral.csv") %>% filter(mobster_test_type == "single")
df_outline <- read_csv("../data/fig5_result_sim_center_edge_neutral_0.25.csv") %>% dplyr::select(sample, seed, location)
df_summary <- df %>% left_join(df_outline, by=c("sample"="sample", "seed"="seed")) %>% 
  group_by(seed, location) %>% summarise(selected_prop = mean(mobster_result_vec == "selected"),.groups = "drop") %>% 
  mutate(outline=ifelse(location == "edge", "edge", "center"))
p <- ggplot(df_summary, aes(x = outline, y = selected_prop, fill=outline)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Proportion of \nselected region detected by MOBSTER", title="Single sample (neutral)") +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif",
                     comparisons = list(c("center", "edge")),
                     size = 7 / .pt)+
  scale_fill_manual(values = c("edge" = "#E69F00", "center" = "#56B4E9")) +
  scale_x_discrete(labels=c("edge"="Edge","center"="Center" ))+
  ylim(0,0.7)+
  my_theme()+
  theme(
    legend.position = "none"
  )
ggsave("figure5h.pdf",unit="cm", width=5, height=5,dpi = 300)


#figure5i----
#Proportion of cases tested positive for natural selection when simulating neutrally evolving tumors.
#When we test using a slide window approach, we observed a higher proportion of selected tests.
df <- read_csv("../data/fig5_mob_single_window_neutral.csv")
df_summary <- df %>%
  group_by(seed, mobster_test_type) %>%
  summarise(selected_prop = mean(mobster_result_vec == "selected"),
            .groups = "drop")
df_paired <- df_summary %>%
  group_by(seed) %>%
  filter(all(c("single","window") %in% mobster_test_type)) %>%
  ungroup()
max_y <- max(df_paired$selected_prop, na.rm = TRUE)
p <- ggplot(df_paired, aes(x = mobster_test_type, y = selected_prop, fill = mobster_test_type)) +
  geom_point(aes(color= mobster_test_type), size=0.5)+
  geom_line(aes(group = seed), colour = "black", linewidth = 0.2, alpha = 0.8) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_manual(values = c("window" = "#E69F00", "single" = "#56B4E9")) +
  scale_fill_manual(values = c("window" = "#E69F00", "single" = "#56B4E9")) +
  scale_x_discrete(
    labels = c("single" = "Single\nsample", "window" = "Sliding\nwindow")
  ) +
  labs(
    x = NULL,
    y = "Proportion of selected results",
    title = "Neutral"
  ) +
  my_theme() +
  theme(
    legend.position = "none"
  ) +
  ylim(0, 0.8)+
  stat_compare_means(
    comparisons = list(c("single","window")),
    method      = "wilcox.test",
    paired      = TRUE,
    label       = "p.signif",
    label.y     = max_y + 0.05
  )
ggsave("figure5i.pdf", unit="cm", width=5, height=5,dpi = 300)


#figure5l----
#Proportion of sectors that are tested positive for natural selection for samples near block boundary and within block for HCC and LUAD
df <- read_csv("../data/fig5_mob_boundary_real_data.csv")
df1 <- df %>%
  mutate(
    boundary_flag = case_when(
      is.logical(boundary) ~ boundary,
      is.numeric(boundary) ~ boundary != 0,
      TRUE ~ tolower(as.character(boundary)) %in% c("true","t","1","yes","y")
    ),
    selected_flag = mobster_result == "selected")

df_summary <- df1 %>%
  group_by(patient, input_type) %>%
  summarise(
    n_boundary              = sum(boundary_flag, na.rm = TRUE),
    n_nonboundary           = sum(!boundary_flag, na.rm = TRUE),
    boundary_selected       = sum(selected_flag &  boundary_flag, na.rm = TRUE),
    nonboundary_selected    = sum(selected_flag & !boundary_flag, na.rm = TRUE),
    boundary_selected_ratio    = ifelse(n_boundary > 0,    boundary_selected    / n_boundary,    NA_real_),
    nonboundary_selected_ratio = ifelse(n_nonboundary > 0, nonboundary_selected / n_nonboundary, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    tumor_type = ifelse(grepl("DT", patient),
                        "liver_cancer", "lung_cancer"))
df_long <- df_summary %>%
  pivot_longer(
    c(boundary_selected_ratio, nonboundary_selected_ratio),
    names_to  = "boundary_type",
    values_to = "selected_ratio"
  ) %>%
  mutate(
    boundary = factor(ifelse(boundary_type == "boundary_selected_ratio", "Boundary", "Non-boundary"),
                      levels = c("Non-boundary", "Boundary"))
  )
df_long_plot <- df_long %>% mutate(boundary=ifelse( boundary == "Boundary", "block boundary", "block interior"))

p <- ggplot(df_long_plot , aes(x = boundary, y = selected_ratio, fill = boundary)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "Proportion of selected\nregion detected by MOBSTER", title="Single sample") +
  scale_fill_manual(values = c("block boundary" = "#E69F00", "block interior" = "#56B4E9")) +
  scale_x_discrete(labels=c("block boundary" = "Block\nboundary", "block interior" = "Block\ninterior"))+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y = 1,
                     comparisons = list(c("block boundary", "block interior")),
                     size = 6 / .pt) +
  ylim(0,1.1)+
  my_theme()+
  theme(
    legend.position = "none",
    strip.text = element_text(size=7),
    plot.title = element_text(face="plain")
  )
ggsave("figure5l.pdf",unit="cm", width=3.5, height=4,dpi = 300)

##figure5m----
#Proportion of windows that are tested positive for natural selection for composite samples with mixed genetic class (i.e. mixed) 
#and pure genetic class (i.e. pure) for HCC and LUAD respectively.
df <- read_csv("../data/fig5_mob_pileup_result_ccf_vaf_window_single.csv")
df_summary <- df %>%
  mutate(type_combo = paste(sample_type, mobster_result)) %>% filter(mobster_test_type == "window", input_type == "ccf_mean") %>%
  group_by(patient, grid_dist) %>%
  summarise(
    pure_select   = sum(type_combo == "pure selected"),
    mixed_select  = sum(type_combo == "mixed selected"),
    pure_neutral    = sum(type_combo == "pure neutral"),
    mixed_neutral   = sum(type_combo == "mixed neutral"),
    .groups = "drop"
  ) %>% filter(grid_dist != -1 , grid_dist != "grid_dist")

df2 <- df_summary %>%
  mutate(
    mixed_denom = mixed_select + mixed_neutral,
    pure_denom  = pure_select  + pure_neutral,
    mixed_prop  = ifelse(mixed_denom > 0, mixed_select / mixed_denom, NA_real_),
    pure_prop   = ifelse(pure_denom  > 0, pure_select  / pure_denom,  NA_real_)
  ) %>%
  dplyr::select(patient, grid_dist, mixed_prop, pure_prop, ) %>% na.omit() %>%
  mutate(tumor_type = ifelse(grepl("DT", patient), "liver_cancer", "lung_cancer"))

df_long <- df2 %>%
  pivot_longer(cols = c(mixed_prop, pure_prop),
               names_to = "prop_type",
               values_to = "value")

df_long_plot <- df_long %>% filter(grid_dist == 0.75) %>% mutate(prop_type = ifelse(prop_type== "mixed_prop", "mixed region", "pure region"))

p2 <- ggplot(df_long_plot, aes(x = prop_type, y = value, fill = prop_type)) +
  geom_boxplot() +
  stat_compare_means(
    comparisons = list(c("mixed region", "pure region")),
    method = "wilcox.test",
    size = 6 / .pt,
    label = "p.signif"
  ) +
  scale_fill_manual(values = c("mixed region" = "#E69F00", "pure region" = "#56B4E9")) +
  scale_x_discrete(labels=c("mixed region" = "Mixed\nregion", "pure region" = "Pure\nregion"))+
  labs(x = NULL, y = "Proportion of selected \nregion detected by MOBSTER", title="Sliding window") +
  ylim(0,1.15)+
  my_theme()+
  theme(
    legend.position = "none",
    strip.text = element_text(size=7),
    plot.title = element_text(face="plain")
  )
ggsave("figure5m.pdf",unit="cm", width=3.5, height=4,dpi = 300)

###figure5n----
df <- read_csv("../data/fig5_mob_boundary_real_data.csv")
df1 <- df %>%
  mutate(
    boundary_flag = case_when(
      is.logical(boundary) ~ boundary,
      is.numeric(boundary) ~ boundary != 0,
      TRUE ~ tolower(as.character(boundary)) %in% c("true","t","1","yes","y")
    ),
    selected_flag = mobster_result == "selected"
  )
df_center_edge <- read_csv("../data/fig5_result_center_edge_all0.25.csv") %>% dplyr::select(location, sample)
df_all <- df1 %>% left_join(df_center_edge, by=c("Punch"="sample")) %>% dplyr::select(mobster_result, patient, Punch, location)
df_long_plot <- df_long %>% filter(input_type %in% c("ccf_mean")) %>% mutate(boundary=ifelse( boundary == "Boundary", "block boundary", "block interior"))
df_prop <- df_all %>%
  group_by(patient, location) %>%
  summarise(
    selected_prop = mean(mobster_result == "selected"),
    .groups = "drop"
  )
ggplot(df_prop, aes(x = location, y = selected_prop, fill = location)) +
  geom_boxplot(outlier.shape = NA) +
  stat_compare_means(
    comparisons = list(c("center", "edge")),
    method = "wilcox.test",
    label = "p.signif",
    size = 6 / .pt
  ) +
  scale_fill_manual(values = c("center" = "#E69F00", "edge" = "#56B4E9")) +
  scale_x_discrete(labels=c("center" = "Center", "edge" = "Edge"))+
  labs(
    title = "Single sample",
    x = "",
    y = "Proportion of selected samples"
  ) +
  ylim(0,1.1)+
  my_theme() +
  theme(
    legend.position = "none",
    plot.margin = margin(1,0,-5,0)
  )
ggsave("figure5n.pdf",unit="cm", width=3, height=4,dpi = 300)

##figure5o----
#dN/dS values calculated for sectors testing positive for natural selection (i.e. selected sample) as well as testing negative
df_dnds <- read_csv("../data/fig5_result_dnds_singles_sample.csv")
df_mobster <- read_csv("../data/fig5_mob_pileup_result_ccf_vaf_window_single.csv")
df_mobster_single <- df_mobster %>% filter(sample_type == "one_sample", input_type=="ccf_mean") %>%
  mutate(samples = gsub("LUAD", "", samples))
df_dnds_all <- df_dnds %>% filter(name == "wall") %>% na.omit()
df_mobster_dnds <- df_dnds_all %>% left_join(df_mobster_single, by=c("SampleID" = "samples")) %>% na.omit() 
comparisons_list <- list(c("selected", "neutral"))
ggplot(df_mobster_dnds, aes(x=mobster_result, y=mle, fill=mobster_result)) +
  geom_boxplot(outlier.shape =NA) +
  stat_compare_means(
    comparisons = comparisons_list,
    method = "wilcox.test",
    paired = FALSE,
    label = "p.signif",
    size = 6 / .pt
  ) + scale_fill_manual(values = c("selected" = "#E69F00", "neutral" = "#56B4E9")) +
  scale_x_discrete(labels=c("selected" = "Selected\nsample", "neutral"="Neutral\nsample"))+
  labs(
    x = NULL,
    y = expression( italic(d["N"]/d["S"]) )
  )+
  ylim(0.5,1.7)+
  my_theme()+
  theme(
    legend.position = "none")
ggsave("figure5p.pdf",unit="cm", width=3.5, height=4,dpi = 300)
