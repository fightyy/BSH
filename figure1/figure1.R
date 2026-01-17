library(tidyverse)
library(this.path)
library(ggpubr)
library(ggpmisc)
library(car)
setwd(this.dir())
source("~/yy_application/R/bin/theme_setup_simple.R")

# read slope data 
df <- read_csv("../data/fig1_result_slope_all_relative.csv")%>% mutate(
  Patient = map_chr(str_split(Patient, "_"), dplyr::last)) %>% filter(Tumor_type != "esophageal_cancer")


###Figure1d----
#The linear relationship between relative (physical) distance and genetic divergence in DT06
fit <- lm(ITH_values ~ Relative_physical_distance, data = df %>% filter(Patient == "DT06"))
r2 <- summary(fit)$r.squared
pval <- summary(fit)$coefficients["Relative_physical_distance", "Pr(>|t|)"]

#Adopt scientific notation
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
r2_str   <- format(round(r2, 2), nsmall = 3)
text_label_expr <- bquote(
  "p" == .(mant_str) %*% 10^.(dp$expo) * "," ~ R^2 == .(r2_str)
)
text_label_str <- as.character(as.expression(text_label_expr))

#plot
liver_plot <- ggplot(df %>% filter(Patient == "DT06"), aes(x = Relative_physical_distance, y = ITH_values)) +
  geom_point(aes(color = Tumor_type), size = 1, alpha = 0.5) +
  geom_smooth(stat="smooth", aes(group = interaction(Patient, Tumor_type), color = Tumor_type), method = "lm", formula= y~x, se = TRUE, linewidth = 0.8, color = "red", linetype = 2) +
  annotate("text", x = mean(range(df %>% filter(Patient == "DT06") %>% .$Relative_physical_distance)), y = max(df$ITH_values),
           label = text_label_str, parse = TRUE,
           size = 7/.pt, fontface = "plain", hjust = 0.5) +
  scale_color_manual(
    name = "Cancer Type",  
    values = c("lung_cancer" = "#1f77b4", "liver_cancer" = "#ff7f0e", "colorectal_cancer" = "#2ca02c", "gastroesophageal_cancer"="#E64B35FF"),
    labels=c("lung_cancer"= "LUAD" ,"liver_cancer"="HCC","colorectal_cancer"="CRC","gastroesophageal_cancer"="GOA")
  ) +
  labs(
    title = "DT06 (HCC)",
    x = "Relative distance",
    y = "Genetic divergence",
    color = "Tumor Type"
  )+
  my_theme()+
  guides(
    color = guide_legend(title.position = "top", title.hjust = 0.5)  
  )+
  theme(
    legend.position = "none"
  )

ggsave("figure1d.pdf",unit="cm", width=4.5, height=5,dpi = 300)


###Figure1e----
#Linear relationships between relative (physical) distance and genetic divergence across all patients
#rename tumor type
df_plot_slope <- df %>% mutate(
  Patient = paste0(Tumor_type, "_", Patient ),
  Tumor_type = dplyr::recode(
    Tumor_type,
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
  Tumor_type = factor(
    Tumor_type,
    levels = c("HCC", "LUAD", "CRC", "GOA", "ESCC", "BLCA", "PCa", "Panc-NEC",  "GC")
  )
)
#plot
ggplot(df_plot_slope, aes(x = Relative_physical_distance, y = ITH_values)) +
  geom_line(stat="smooth", aes(group = interaction(Patient, Tumor_type), color = Tumor_type), method = "lm", formula= y~x,  se = FALSE, linewidth = 0.5,  alpha=0.5, fullrange = TRUE) +
  expand_limits(x = 0)+
  coord_cartesian(xlim = c(0, max(df_plot_slope$Relative_physical_distance)))+
  scale_color_manual(
    name = "Tumor type ",
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
  ) +
  labs(
    x = "Relative distance",
    y = "Genetic divergence",
    color = "Tumor type"
  )+
  my_theme()+
  guides(
    color = guide_legend(title.position = "top", title.hjust = 0.5)  
  )+
  theme(
    legend.position = "right",
    legend.box.margin = margin(t = 0, unit = "cm"),
    plot.margin = margin(t = 0.1, unit = "cm")
  )
ggsave("figure1e.pdf",unit="cm", width=8, height=5,dpi = 300)


###figure 1f----
#Boxplot of regression slopes of the linear relationship between relative physical distance and genetic divergence
#rename tumor type
df_boxplot <- read_csv("../data/fig1_result_slope_all_relative_summary.csv") %>% filter(tumor_type != "esophageal_cancer") %>% 
  mutate(patient = map_chr(str_split(patient, "_"), dplyr::last)) %>% mutate(
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
      #levels = c("LUAD","GOA","GC","EC","CRC","PC","BC","PrC","HCC"),
      levels = c("HCC", "LUAD", "CRC", "GOA", "BLCA", "PCa", "Panc-NEC",  "GC")
    )
  )

#Compare slope pairwise and calculate the p-value
comparisons_list <- combn(levels(df_boxplot$tumor_type), 2, simplify = FALSE)
p_dat <- compare_means(
  slope ~ tumor_type,
  data   = df_boxplot,
  method = "wilcox.test",
  paired = FALSE
)

#transfer p value into * 
p_dat <- p_dat %>%
  mutate(label = cut(p,
                     breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                     labels = c("****","***","**","*","ns")))

# 3) order comparisons_list 
order_g1 <- vapply(comparisons_list, `[[`, "", 1)
order_g2 <- vapply(comparisons_list, `[[`, "", 2)

p_dat <- p_dat %>%
  filter(group1 %in% order_g1 & group2 %in% order_g2) %>%
  arrange(match(paste(group1, group2), paste(order_g1, order_g2)))

# Assign the y position to each comparison in sequence
#ymax <- max(df_boxplot$slope, na.rm = TRUE)
ymax <- 0.7
ymin <- min(df_boxplot$slope, na.rm = TRUE)
step <- 0.1 * (ymax - ymin)
p_dat$y.position <- ymax + step * seq_len(nrow(p_dat))

#plot
p <- ggplot(df_boxplot, aes(x = tumor_type, y = slope, fill = tumor_type)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(
    name = "Tumor type ",
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
  ) +
  labs(title = NULL,
       x = NULL,
       y = "Slope") +
  my_theme() +
  theme(
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_blank()) 
ggsave("figure1f.pdf", unit="cm", width=6, height=5, dpi=300)

###figure1g----
#Boxplot of the R2 statistic from the linear regression for patients with negative slope (Slope <0) and positive slope (Slope <0)
df_boxplot$slope_group <- ifelse(df_boxplot$slope > 0, "Slope > 0", "Slope < 0")
ggplot(df_boxplot, aes(x = slope_group, y = r_squared, fill = slope_group)) +
  geom_boxplot(outlier.shape = 21) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     label.y = max(df_boxplot$r_squared) * 1.05,
                     comparisons = list(c("Slope > 0", "Slope < 0"))) +
  scale_fill_manual(values=c( "Slope > 0"       = "#E64B35FF",  
                             "Slope < 0"      = "#1f77b4" ))+
  labs(x = "Slope group", y = expression(R^2)) +
  ylim(0,1.1)+
  my_theme()+
  theme(legend.position = "none")
ggsave("figure1g.pdf",unit="cm", width=4, height=5,dpi = 300)


###figure1h----
#Boxplot of the number of sampled sectors for patients with negative slope (Slope <0) and positive slope (Slope <0)
df_boxplot <- df_boxplot %>%
  mutate(real_sample_num = (1 + sqrt(1 + 8 * sample_number)) / 2)
ggplot(df_boxplot, aes(x = slope_group, y = real_sample_num, fill = slope_group)) +
  geom_boxplot(outlier.shape = 21) +
  #geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     label.y = max(df_boxplot$real_sample_num) * 1.05,
                     comparisons = list(c("Slope > 0", "Slope < 0"))) +
  scale_fill_manual(values=c( "Slope > 0"       = "#E64B35FF",  
                              "Slope < 0"      = "#1f77b4" ))+
  labs(x = "Slope group", y = "Sample number") +
  ylim(0,35)+
  my_theme()+
  theme(legend.position = "none")
ggsave("figure1h.pdf",unit="cm", width=4, height=5,dpi = 300)


###figure1i----
#Proportion of total variation in the regression slopes explained by variance between different patients or between tumors
fit <- lm(slope ~ tumor_type, data = df_boxplot)
anova_res <- car::Anova(fit, type = 2)

between_var <- anova_res["tumor_type", "Sum Sq"]
within_var  <- anova_res["Residuals", "Sum Sq"]

total_var   <- between_var + within_var
between_prop <- between_var / total_var
within_prop  <- within_var / total_var

var_df <- data.frame(
  Source   = c("Between tumor type", "Within tumor type"),
  Variance = c(between_prop, within_prop)
)

ggplot(var_df, aes(x = Source, y = Variance, fill = Source)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c(
    "Between tumor type" = "#E64B35FF",
    "Within tumor type"  = "#4DBBD5FF"
  )) +
  scale_x_discrete(labels = c(
    "Between tumor type" = "Between\ntumor type",
    "Within tumor type"  = "Within\ntumor type"
  )) +
  labs(
    y = "Proportion of slope variance explained",
    x = NULL
  ) +
  my_theme() +
  theme(
    legend.position = "none"
  )

ggsave("figure1i.pdf",unit="cm", width=4, height=5,dpi = 300)
 

