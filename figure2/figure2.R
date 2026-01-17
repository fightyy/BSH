library(tidyverse)
library(this.path)
library(adegenet)  # chooseCN, monmonier
library(spdep)     # nb2listw/nb2mat（可选）
library(graphics)  # plot
library(ape)
library(phangorn)
library(ggpmisc)
library(patchwork)
library(spdep)
library(adegenet)
library(ggtree)
setwd(this.dir())
source("../bin/theme_setup_simple.R")

# get functions needed
source("../bin/function.R")

###figure2a----
#Spatial distribution of tumor boundaries in DT10.
#Circles represent tumor samples, lines connect neighboring samples, and arrows indicate a boundary.
#The background heatmap represents the levels regional ITH along the 2D space.
patient_cancer <- "liver_cancer_DT10"
relative_thresh<- 0.8
#read boundary data
boundary_df <- read_csv(paste0("../data/fig2_monmonier_out_two_path/", patient_cancer, "_", relative_thresh, "_monmonier_result.csv"))

# Color palette used for tree and spatial plots
tree_color <- c("#F8766D","#00BFC4","#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))

# Parse tumor type and patient ID from `patient_cancer`
# e.g. "liver_cancer_DT10" -> tumor_type = "liver_cancer", patient = "DT10"
patient_tumor_vec <- strsplit( patient_cancer, "_")
tumor_type <- paste0(patient_tumor_vec[[1]][1],"_cancer")
patient <- patient_tumor_vec[[1]][3]

# Read spatial coordinates of bulk samples (absolute, non-size-normalised)
location <- read_csv("../data/fig2_public_liver_lung_loaction_absolute_nsr.csv") %>% filter(Patient==patient, Tumor_type==tumor_type) %>%
  dplyr::select(Sample, X, Y ) %>% dplyr::rename(x=X, y=Y, sample=Sample)  %>% mutate(z=0)
r <- max(max(location$x)-min(location$x), max(location$y)-min(location$y))/2

# Read VAF matrix for this patient
mutation_code <- read_csv(paste0("../data/vaf_public_liver_lung/", patient_cancer, "_vaf.csv"))
mutation_code <- mutation_code %>% dplyr::select(-"mut_id")
# Use only samples that appear in both mutation matrix and location table
sample_names <- intersect(colnames(mutation_code), location$sample)
sample_names <- sample_names[sample_names != "NU"]
mutation_code <- mutation_code[,sample_names]
# Add GL (germline) column initialized to 0
mutation_code$GL <- 0
mutation_code = mutation_code > 0
storage.mode(mutation_code) = "numeric"
phy_data=mutation_code %>% t() %>% phangorn::phyDat(type="USER", levels=c(0,1))
#get tree file
rpar = readRDS(paste0("../data/rds_public_liver_lung/",patient_cancer, "_tree.rds"))

# Get best clustering based on CH index
genetic_cluster <- get_CH_index(rpar)
print("###### tree plot finished ######")

# Map cluster labels to spatial coordinates of bulk samples
loc_genetic_cluster <- data.frame()
for(sample_ in sample_names){
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

# Compute continuous ITH surface across tumor
select_info <- loc_genetic_cluster[c("X","Y")]
ITH_matrix <- get_ITH_matrix(phy_data,select_info,r)
print("###### ITH matirx finished ######")

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

#plot boundary
loc <- location %>% mutate(sample = as.character(sample)) %>% distinct(sample, .keep_all = TRUE) %>% filter(sample %in% sample_names)
#get gabriel graph data
nb.gab <- build_nb_gabriel(loc)
xy <- as.matrix(loc[, c("x","y")])
#Convert the adjacent nb into edge data used for drawing
edges <- nb_to_edges_df(nb.gab, xy, samples=loc$sample)
#Convert the Monmonier's  data into boudary data used for drawing
arrow_df <- make_arrow_segments_all(boundary_df)
arrow_df <- arrow_df %>%
  dplyr::rename(Path = run_id) %>%
  mutate(Path = gsub("^run", "Path", Path))

#plot
p1 <- ggplot(all_point_ITH, aes(x = x, y = y)) +
  stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
  labs(fill="ITH")+
  scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                       breaks=c(0,0.25,0.5,0.75,1),colours =  colormap,
                       guide = guide_colorbar(order = 2))+
  geom_segment(data = edges,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.4, alpha = 0.4, colour = "black", lineend = "round") +
  geom_segment(data = arrow_df,
               aes(x = bound_x, y = bound_y, xend = xend, yend = yend,
                   group = interaction(Path, dir), colour = Path),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               linewidth = 0.6, lineend = "round") +
  scale_colour_manual(
    values = c(
      "Path1" = "black",
      "Path2" = "#4393C3"
    ),
    name = "Boundary",
    labels = c(
      "Path1" = "Run1",
      "Path2" = "Run2"
    )
  )+
  labs(x = "X", y = "Y", fill = "ITH")

ITH_plot <- p1 +
  ggnewscale::new_scale_fill()+
  geom_point(data = loc_genetic_cluster,aes(x=X,y=Y), fill = "grey80", size=3.5,shape=21,stroke=0.4)+
  scale_fill_manual(values=tree_color)+labs(fill="class")+
  labs(title = paste0("DT10 (HCC)"), x = NULL, y = NULL, fill="Class") +
  guides(fill = guide_legend(order=1))+
  #ylim(1.4, 4.5)+
  my_theme()+
  theme(
    legend.box.margin = margin(t = -0.2, unit = "cm"),
    plot.margin = margin(t = 0.1, unit = 'cm'))

ggsave("figure1a.pdf", width=6, height=5, unit="cm", dpi = 300 )



##figure2b----
#Spatial distribution of tumor boundaries in DT13

#code comment see figure2a
patient_cancer <- "liver_cancer_DT13"
relative_thresh<- 0.5
boundary_df <- read_csv(paste0("../data/fig2_monmonier_out_two_path/", patient_cancer, "_", relative_thresh, "_monmonier_result.csv"))
tree_color <- c("#F8766D","#00BFC4","#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
patient_tumor_vec <- strsplit( patient_cancer, "_")
tumor_type <- paste0(patient_tumor_vec[[1]][1],"_cancer")
patient <- patient_tumor_vec[[1]][3]
location <- read_csv("../data/fig2_public_liver_lung_loaction_absolute_nsr.csv") %>% filter(Patient==patient, Tumor_type==tumor_type) %>%
  dplyr::select(Sample, X, Y ) %>% dplyr::rename(x=X, y=Y, sample=Sample)  %>% mutate(z=0)
r <- max(max(location$x)-min(location$x), max(location$y)-min(location$y))/2
mutation_code <- read_csv(paste0("../data/vaf_public_liver_lung/", patient_cancer, "_vaf.csv"))
mutation_code <- mutation_code %>% dplyr::select(-"mut_id")
sample_names <- intersect(colnames(mutation_code), location$sample)
sample_names <- sample_names[sample_names != "NU"]
mutation_code <- mutation_code[,sample_names]
mutation_code$GL <- 0
mutation_code = mutation_code > 0
storage.mode(mutation_code) = "numeric"
phy_data=mutation_code %>% t() %>% phangorn::phyDat(type="USER", levels=c(0,1))
rpar = readRDS(paste0("../data/rds_public_liver_lung/",patient_cancer, "_tree.rds"))
genetic_cluster <- get_CH_index(rpar)
loc_genetic_cluster <- data.frame()
for(sample_ in sample_names){
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
print("###### ITH matirx finished ######")

set.seed(123)
ITH <- ITH_matrix$ITH
temp_order <- c()
for (n in 1:length(ITH)){
  if (sample(c(0,1),1,prob=c(1-ITH[n],ITH[n]))==1){
    temp_order=c(temp_order,n)
  }
}
all_point_ITH <- ITH_matrix[temp_order,]

loc <- location %>% mutate(sample = as.character(sample)) %>% distinct(sample, .keep_all = TRUE) %>% filter(sample %in% sample_names)
nb.gab <- build_nb_gabriel(loc)
xy <- as.matrix(loc[, c("x","y")])
edges <- nb_to_edges_df(nb.gab, xy, samples=loc$sample)
arrow_df <- make_arrow_segments_all(boundary_df)
arrow_df <- arrow_df %>%
  dplyr::rename(Path = run_id) %>%
  mutate(Path = gsub("^run", "Path", Path))

p1 <- ggplot(all_point_ITH, aes(x = x, y = y)) +
  stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
  labs(fill="ITH")+
  scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                       breaks=c(0,0.25,0.5,0.75,1),colours =  colormap,
                       guide = guide_colorbar(order = 2))+
  geom_segment(data = edges,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.4, alpha = 0.4, colour = "black", lineend = "round") +
  geom_segment(data = arrow_df,
               aes(x = bound_x, y = bound_y, xend = xend, yend = yend,
                   group = interaction(Path, dir), colour = Path),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               linewidth = 0.6, lineend = "round") +
  scale_colour_manual(
    values = c(
      "Path1" = "black",
      "Path2" = "#4393C3"
    ),
    name = "Boundary",
    labels = c(
      "Path1" = "Run1",
      "Path2" = "Run2"
    )
  )+
  labs(x = "X", y = "Y", fill = "ITH")

ITH_plot <- p1 +
  ggnewscale::new_scale_fill()+
  geom_point(data = loc_genetic_cluster,aes(x=X,y=Y), fill = "grey80", size=3.5,shape=21,stroke=0.4)+
  scale_fill_manual(values=tree_color)+labs(fill="class")+
  labs(title = paste0("DT13 (HCC)"), x = NULL, y = NULL, fill="Class") +
  guides(fill = guide_legend(order=1))+
  my_theme()+
  theme(
    legend.box.margin = margin(t = -0.2, unit = "cm"),
    plot.margin = margin(t = 0.1, unit = 'cm'))


ggsave("figure2b.pdf", width=6, height=5, unit="cm", dpi = 300 )


###figure2c-2d----
#figure2c:Maximum parsimony tree constructed from somatic mutations discovered in sectors of DT10
#figure2d:Two-dimensional spatial distribution of tumor sectors in DT10
#code comment see figure2a
patient_cancer <- "liver_cancer_DT10"
tree_color <- c("#F8766D","#00BFC4","#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
patient_tumor_vec <- strsplit( patient_cancer, "_")
tumor_type <- paste0(patient_tumor_vec[[1]][1],"_cancer")
patient <- patient_tumor_vec[[1]][3]
location <- read_csv("../data/fig2_public_liver_lung_loaction_absolute_nsr.csv") %>% filter(Patient==patient, Tumor_type==tumor_type) %>%
  dplyr::select(Sample, X, Y ) %>% dplyr::rename(x=X, y=Y, sample=Sample)  %>% mutate(z=0)
r <- max(max(location$x)-min(location$x), max(location$y)-min(location$y))/2
mutation_code <- read_csv(paste0("../data/vaf_public_liver_lung/", patient_cancer, "_vaf.csv"))
mutation_code <- mutation_code %>% dplyr::select(-"mut_id")
sample_names <- intersect(colnames(mutation_code), location$sample)
sample_names <- sample_names[sample_names != "NU"]
mutation_code <- mutation_code[,sample_names]
mutation_code$GL <- 0
mutation_code = mutation_code > 0
storage.mode(mutation_code) = "numeric"
phy_data=mutation_code %>% t() %>% phangorn::phyDat(type="USER", levels=c(0,1))
rpar = readRDS(paste0("../data/rds_public_liver_lung/",patient_cancer, "_tree.rds"))
genetic_cluster <- get_CH_index(rpar)

loc_genetic_cluster <- data.frame()
for(sample_ in sample_names){
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
print("###### ITH matirx finished ######")

set.seed(123)
ITH <- ITH_matrix$ITH
temp_order <- c()
for (n in 1:length(ITH)){
  if (sample(c(0,1),1,prob=c(1-ITH[n],ITH[n]))==1){
    temp_order=c(temp_order,n)
  }
}
all_point_ITH <- ITH_matrix[temp_order,]
p1 <- ggplot(all_point_ITH, aes(x = x, y = y)) +
  stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
  labs(fill="ITH")+
  scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                       breaks=c(0,0.25,0.5,0.75,1),colours = colormap)

ITH_plot <- p1 +
  ggnewscale::new_scale_fill()+
  geom_point(
    data = loc_genetic_cluster,
    aes(x = X, y = Y,fill= cluster),
    size=3.5,stroke=0.4,shape=21
  ) +
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=gsub(paste0(patient_tumor_vec[[1]][3],"_"),"", rownames(loc_genetic_cluster))),size=6/.pt,vjust=-1)+
  scale_color_manual(name = "Genetict class", values = tree_color) +
  guides(fill = guide_legend(order=1))+
  labs(title = "DT10 (HCC)", x = NULL, y = NULL, fill="Class") +
  my_theme()+
  theme(
    legend.box.margin = margin(t = -0.2, unit = "cm"),
    plot.margin = margin(t = 0.1, unit = 'cm'))
ggsave("figure2c.pdf", width=6, height=5, unit="cm", dpi = 300 )

group <- genetic_cluster$Genetic_group
# Build annotation data frame mapping tip labels to cluster and sample ID
annotation <- data.frame(label=names(group),group=as.character(group),sample=names(group))
rpar1 <- drop.tip(rpar, "GL")
trs <- full_join(rpar1,annotation,by='label')
xmax <- max(phytools::nodeHeights(rpar1))
scale_width <- 30
p1 <- ggtree(trs) +
  geom_tippoint(aes(color = group), shape = 16, size = 10/.pt) +
  geom_text2(aes(label=gsub("DT10_", "", sample)),size=4/.pt,color="white")+
  geom_treESCCle(x = xmax - scale_width, y = 3, width = scale_width, fontsize = 6/.pt)+
  scale_color_manual(values = tree_color, name="Genetic class")+
  geom_rootedge(rootedge = rpar$edge.length[length(rpar$edge.length)]) +
  coord_cartesian(clip = 'off') +
  geom_rootpoint()+
  labs(title="DT10 (HCC)")+
  my_theme()+
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
ggsave("figure2d.pdf", unit="cm", width=5, height=5,dpi = 300)



###figure2e-2f
# figure2e:Maximum parsimony tree constructed from somatic mutations discovered in sectors of DT13
# figure2f:Two-dimensional spatial distribution of tumor sectors in DT13
patient_cancer <- "liver_cancer_DT13"
tree_color <- c("#F8766D","#00BFC4","#C77CFF","#7CAE00","#3288BD", "#FFFFBF","#FDAE61", "#9E0142","#5E4FA2")
colormap=rev(RColorBrewer::brewer.pal(11,'Spectral'))
patient_tumor_vec <- strsplit( patient_cancer, "_")
tumor_type <- paste0(patient_tumor_vec[[1]][1],"_cancer")
patient <- patient_tumor_vec[[1]][3]
location <- read_csv("../data/fig2_public_liver_lung_loaction_absolute_nsr.csv") %>% filter(Patient==patient, Tumor_type==tumor_type) %>%
  dplyr::select(Sample, X, Y ) %>% dplyr::rename(x=X, y=Y, sample=Sample)  %>% mutate(z=0)
r <- max(max(location$x)-min(location$x), max(location$y)-min(location$y))/2
mutation_code <- read_csv(paste0("../data/vaf_public_liver_lung/", patient_cancer, "_vaf.csv"))
mutation_code <- mutation_code %>% dplyr::select(-"mut_id")
sample_names <- intersect(colnames(mutation_code), location$sample)
sample_names <- sample_names[sample_names != "NU"]
mutation_code <- mutation_code[,sample_names]
mutation_code$GL <- 0
mutation_code = mutation_code > 0
storage.mode(mutation_code) = "numeric"
phy_data=mutation_code %>% t() %>% phangorn::phyDat(type="USER", levels=c(0,1))
rpar = readRDS(paste0("../data/rds_public_liver_lung/",patient_cancer, "_tree.rds"))
genetic_cluster <- get_CH_index(rpar)

loc_genetic_cluster <- data.frame()
for(sample_ in sample_names){
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
print("###### ITH matirx finished ######")

set.seed(123)
ITH <- ITH_matrix$ITH
temp_order <- c()
for (n in 1:length(ITH)){
  if (sample(c(0,1),1,prob=c(1-ITH[n],ITH[n]))==1){
    temp_order=c(temp_order,n)
  }
}
all_point_ITH <- ITH_matrix[temp_order,]
p1 <- ggplot(all_point_ITH, aes(x = x, y = y)) +
  stat_density_2d(geom = "raster",aes(fill = ..ndensity..),contour = F)+
  labs(fill="ITH")+
  scale_fill_gradientn(labels=round(seq(0,max(all_point_ITH$ITH)/1,length.out=5),digits = 2),
                       breaks=c(0,0.25,0.5,0.75,1),colours = colormap)

ITH_plot <- p1 +
  ggnewscale::new_scale_fill()+
  geom_point(
    data = loc_genetic_cluster,
    aes(x = X, y = Y,fill= cluster),
    size=3.5,stroke=0.4,shape=21
  ) +
  geom_text(data = loc_genetic_cluster,aes(x=X,y=Y,label=gsub(paste0(patient_tumor_vec[[1]][3],"_"),"", rownames(loc_genetic_cluster))),size=6/.pt,vjust=-1)+
  scale_color_manual(name = "Genetict class", values = tree_color) +
  guides(fill = guide_legend(order=1))+
  labs(title = "DT13 (HCC)", x = NULL, y = NULL, fill="Class") +
  my_theme()+
  theme(
    legend.box.margin = margin(t = -0.2, unit = "cm"),
    plot.margin = margin(t = 0.1, unit = 'cm'))
ggsave("figure2e.pdf", width=6, height=5, unit="cm", dpi = 300 )

group <- genetic_cluster$Genetic_group
annotation <- data.frame(label=names(group),group=as.character(group),sample=names(group))
rpar1 <- drop.tip(rpar, "GL")
trs <- full_join(rpar1,annotation,by='label')
xmax <- max(phytools::nodeHeights(rpar1))
scale_width <- 30
p1 <- ggtree(trs) +
  geom_tippoint(aes(color = group), shape = 16, size = 10/.pt) +
  geom_text2(aes(label=gsub("DT13_", "", sample)),size=4/.pt,color="white")+
  geom_treESCCle(x = xmax - scale_width, y = 3, width = scale_width, fontsize = 6/.pt)+
  scale_color_manual(values = tree_color, name="Genetic class")+
  geom_rootedge(rootedge = rpar$edge.length[length(rpar$edge.length)]) +
  coord_cartesian(clip = 'off') +
  geom_rootpoint()+
  labs(title="DT13 (HCC)")+
  my_theme()+
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
p1
ggsave("figure2f.pdf", unit="cm", width=5, height=5,dpi = 300)

###figure2h----
#Boxplot of ENNMN values across tumor types
df <- read_csv("../data/fig2_result_enn_treestat_real_data.csv")
cancer_type <- unique(df$tumor_type)
df_boxplot <- df %>%
  mutate(
    tumor_type = dplyr::recode(
      tumor_type,
      lung_cancer             = "LUAD",
      gastroesophageal_cancer = "GOA",
      gastric_cancer          = "GC",
      esophageal_cancer       = "ESCC",
      colorectal_cancer       = "CRC",
      pancreatic_cancer       = "Panc-NEC",
      bladder_cancer          = "BLCA",
      prostate_cancer         = "PCa",
      liver_cancer            = "HCC"
    ),
    tumor_type = factor(
      tumor_type,
      levels = c("HCC", "LUAD", "CRC", "GOA", "ESCC", "BLCA", "PCa", "Panc-NEC",  "GC")
    )
  )
comparisons_list <- combn(unique(df_boxplot$tumor_type), 2, simplify = F)
ggplot(df_boxplot, aes(x = tumor_type, y = value_enn, fill = tumor_type)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "HCC"       = "#E64B35FF",
      "LUAD"      = "#1f77b4",
      "CRC"       = "#2ca02c",
      "GOA"       = "#ff7f0e",
      "ESCC"      = "#FBB4AE",
      "BLCA"      = "#B3CDE3",
      "PCa"       = "#CCEBC5",
      "Panc-NEC"  = "#DECBE4",
      "GC"        = "#FED9A6"
    )
  )+

  labs(
    x = NULL,
    y = expression(ENN[MN])) +
  my_theme()+
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank(),
  )
ggsave("figure2h.pdf",unit="cm", width=7, height=5,dpi = 300)

#figure2i----
#Distribution of the Silhouette values across all patients
df <- read_csv("../data/fig2_result_boundary_real.csv")
ggplot(df, aes(x = silhouette_mean)) +
  geom_histogram(
    aes(y = after_stat(count)),
    binwidth = 0.1,
    boundary = 0,
    fill = "#1f77b4", color = "black"
  ) +
  labs(
    x = "Silhouette value",
    y = "Count"
  ) +
  my_theme()+
  theme(
    axis.text.x = element_text(size = 6)
  )
ggsave("figure2i.pdf", unit="cm", width=6, height=5, dpi=300)

###figure2j----
#Boxplot of Silhouette values across tumor types
df <- read_csv("../data/fig2_result_boundary_real.csv")
df_boxplot <-  df %>% 
  mutate(patient = map_chr(str_split(Patient, "_"), dplyr::last)) %>% mutate(
    tumor_type = dplyr::recode(
      tumor_type,
      lung_cancer             = "LUAD",
      gastroesophageal_cancer = "GOA",
      gastric_cancer          = "GC",
      esophageal_cancer       = "ESCC",
      colorectal_cancer       = "CRC",
      pancreatic_cancer       = "Panc-NEC",
      bladder_cancer          = "BLCA",
      prostate_cancer         = "PCa",
      liver_cancer            = "HCC"
    ),
    # 如果你还想保证画图时按照这个顺序出现：
    tumor_type = factor(
      tumor_type,
      levels = c("HCC", "LUAD", "CRC", "ESCC", "GOA", "BLCA", "PCa", "Panc-NEC",  "GC")
    )
  )
comparisons_list <- combn(levels(df_boxplot$tumor_type), 2, simplify = FALSE)

p_dat <- compare_means(
  silhouette_mean ~ tumor_type,
  data   = df_boxplot,
  method = "wilcox.test",
  paired = FALSE
)
p_dat <- p_dat %>%
  mutate(label = cut(
    p,
    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
    labels = c("****","***","**","*","ns")
  ))
# Only retain the comparisons with significant differences
p_dat <- p_dat %>% filter(label != "ns")
order_g1 <- vapply(comparisons_list, `[[`, "", 1)
order_g2 <- vapply(comparisons_list, `[[`, "", 2)

p_dat <- p_dat %>%
  filter(group1 %in% order_g1 & group2 %in% order_g2) %>%
  arrange(match(paste(group1, group2), paste(order_g1, order_g2)))

ymax <- 0.85
ymin <- min(df_boxplot$silhouette_mean, na.rm = TRUE)
step <- 0.1 * (ymax - ymin)
p_dat$y.position <- ymax + step * seq_len(nrow(p_dat))

p <- ggplot(df_boxplot, aes(x = tumor_type, y = silhouette_mean, fill = tumor_type)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(
    name = "Tumor type ",
    values = c(
      "HCC"      = "#E64B35FF",
      "LUAD"     = "#1f77b4",
      "CRC"      = "#2ca02c",
      "GOA"      = "#ff7f0e",
      "ESCC"     = "#FBB4AE",
      "BLCA"     = "#B3CDE3",
      "PCa"      = "#CCEBC5",
      "Panc-NEC" = "#DECBE4",
      "GC"       = "#FED9A6"
    )
  ) +
  stat_pvalue_manual(
    data = p_dat,     
    label = "label",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01, bracket.size = 0.3,
    inherit.aes = FALSE, size = 7/.pt
  ) +
  labs(x = NULL, y = "Silhouette value") +
  my_theme() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_blank()
  )
ggsave("figure2j.pdf", unit="cm", width=7, height=5, dpi=300)

