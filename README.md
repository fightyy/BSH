### BSH

This is the code for **Block-Shaped Heterogeneity as an Emergent Principle of Spatial Heterogeneity Across Human Solid Tumors**.

#### *Directory structure*

The code for BSH is organized into different directories and scripts.

The directory structure is as follows:

- bin :  This directory contains ggplot2 theme and function that is needed for regenerating plots.
- data : This directory contains data that is needed for regenerating plots and descriptions for the data .
- figure1-5 : This directory contains the analysis scripts for regenerating main figures  .
- simulation :  This directory contains the deme-based simulation scripts and descriptions for the usage .

#### *Description of Scripts and How to Use*

##### figure1-5

- Description: These folders contains R scripts to generate plots used in the article. All data can be found in the `data` directory. You can reproduce the results by running these R scripts like `Rscript figure1.R`  .

##### simulation

Description: These folders contains python scripts to simulate deme-baed tumor. The simulation parameters and instructions for use can be found in the README file in this directory.

#### *Dependencies*

R session info:

```R
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.6.1

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
  [1] ade4_1.7-22             tidyselect_1.2.0        dplyr_1.1.4             phytools_2.1-1          optimParallel_1.0-2    
  [6] fastmap_1.1.1           combinat_0.0-8          promises_1.2.1          digest_0.6.34           mime_0.12              
 [11] lifecycle_1.0.4         sf_1.0-16               cluster_2.1.6           ellipsis_0.3.2          magrittr_2.0.3         
 [16] compiler_4.3.1          rlang_1.1.3             tools_4.3.1             igraph_2.0.2            utf8_1.2.4             
 [21] phangorn_2.11.1         clusterGeneration_1.3.8 htmlwidgets_1.6.4       pkgbuild_1.4.3          sp_2.1-3               
 [26] classInt_0.4-10         mnormt_2.1.1            scatterplot3d_0.3-44    plyr_1.8.9              RColorBrewer_1.1-3     
 [31] pkgload_1.3.4           KernSmooth_2.23-22      miniUI_0.1.1.1          expm_0.999-9            purrr_1.0.2            
 [36] numDeriv_2016.8-1.1     grid_4.3.1              fansi_1.0.6             urlchecker_1.0.1        profvis_0.3.8          
 [41] xtable_1.8-4            e1071_1.7-14            colorspace_2.1-0        ggplot2_3.5.0           scales_1.3.0           
 [46] iterators_1.0.14        MASS_7.3-60.0.1         cli_3.6.2               vegan_2.6-6.1           remotes_2.4.2.1        
 [51] generics_0.1.3          rstudioapi_0.15.0       reshape2_1.4.4          sessioninfo_1.2.2       spdep_1.3-13           
 [56] cachem_1.0.8            DBI_1.2.2               ape_5.7-1               proxy_0.4-27            stringr_1.5.1          
 [61] splines_4.3.1           maps_3.4.2              parallel_4.3.1          s2_1.1.7                devtools_2.4.5         
 [66] vctrs_0.6.5             boot_1.3-30             Matrix_1.6-5            spData_2.3.4            seqinr_4.2-36          
 [71] foreach_1.5.2           ggnewscale_0.4.10       units_0.8-5             glue_1.7.0              codetools_0.2-19       
 [76] stringi_1.8.3           gtable_0.3.4            later_1.3.2             deldir_2.0-4            quadprog_1.5-8         
 [81] munsell_0.5.0           tibble_3.2.1            pillar_1.9.0            htmltools_0.5.7         R6_2.5.1               
 [86] wk_0.9.2                doParallel_1.0.17       shiny_1.8.0             lattice_0.22-5          memoise_2.0.1          
 [91] httpuv_1.6.14           class_7.3-22            adegenet_2.1.10         Rcpp_1.0.12             fastmatch_1.1-4        
 [96] coda_0.19-4.1           nlme_3.1-164            permute_0.9-7           mgcv_1.9-1              usethis_2.2.3          
[101] fs_1.6.3                pkgconfig_2.0.3
```

Python packages：

```python
## Requirements
- Python **3.x**
- numpy
- sys
- math
- random
- heapq
- subprocess
- collections
```

#### *Data source*

The dataset comprises 165 patients, including private sequencing data generated for this study (five HCC patients, DT42–DT51, and eight LUAD patients). The remaining 152 patients were collected from previously published cohorts, as listed in /data/fig1_patientID_clean.csv. For all 165 patients, VAF matrices are stored in /data/vaf_public_liver_lung, maximum parsimony trees are stored in /data/rds_public_liver_lung, and two-dimensional coordinate files are stored in fig2_public_liver_lung_loaction_absolute_nsr.csv.

