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

- Description: These folders contain R scripts to generate plots used in the article. All data can be found in the `data` directory. You can reproduce the results by running these R scripts like `Rscript figure1.R`  .

##### simulation

Description: This folder contains python scripts to simulate deme-based tumor. The simulation parameters and instructions for use can be found in the README file in this directory.

#### *Dependencies*

R packages：

```R
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.6.1
### R Dependencies
|Package          |Version |
|:----------------|:-------|
|adegenet         |2.1.10  |
|ape              |5.7.1   |
|car              |3.1.2   |
|cowplot          |1.1.3   |
|GET              |1.0.4   |
|ggplot2          |3.5.2   |
|ggpmisc          |0.5.5   |
|ggpubr           |0.6.0   |
|ggtext           |0.1.2   |
|ggtree           |3.10.1  |
|graphics         |4.3.1   |
|ks               |1.14.2  |
|landscapemetrics |2.2.1   |
|mobster          |1.0.0   |
|patchwork        |1.2.0   |
|phangorn         |2.11.1  |
|purrr            |1.0.2   |
|raster           |3.6.26  |
|spatstat         |3.3.0   |
|spdep            |1.3.13  |
|terra            |1.7.78  |
|this.path        |2.4.0   |
|tidyverse        |2.0.0   |
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

The dataset comprises 165 patients, including private sequencing data generated for this study (five HCC patients, DT42–DT51, and eight LUAD patients). The remaining 152 patients were collected from previously published cohorts, as listed in /data/fig1_patientID_clean.csv. For all 165 patients, VAF matrices are stored in /data/vaf_public_liver_lung, maximum parsimony trees are stored in /data/rds_public_liver_lung, and two-dimensional coordinate files are stored in /data/fig2_public_liver_lung_loaction_absolute_nsr.csv.

