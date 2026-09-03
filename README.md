### BSH

This is the code for **block-shaped heterogeneity as an emergent organizing principle of spatial heterogeneity across human solid tumors**.

#### *Directory structure*

The code for BSH is organized into different directories and scripts.

The directory structure is as follows:

- bin :  This directory contains ggplot2 theme and function that is needed for regenerating plots.
- data : This directory contains data that is needed for regenerating plots and descriptions for the data .
- figure1-5 : This directory contains the analysis scripts for regenerating main figures  .
- simulation :  This directory contains the deme-based simulation scripts and descriptions for the usage .

#### *Description of Scripts and How to Use*

##### figure1-5

- Description: These folders contain R scripts to generate main plots in the article. All data can be found in the `data` directory. You can reproduce the results by running these R scripts like `Rscript figure1.R`  .

##### simulation

Description: This folder contains python scripts to simulate deme-based tumor. The simulation parameters, instructions for use and expected output can be found in the README file in this directory.

#### *System requirements and dependencies*

The code was developed and tested on the following environment:

- macOS Sonoma 14.6.1
- R version 4.3.1 
- Platform: aarch64-apple-darwin20 (64-bit)

No non-standard hardware is required. The analyses can be performed on a standard desktop computer.

##### R dependencies

The following R packages are required:

| Package          | Version |
| :--------------- | :------ |
| adegenet         | 2.1.10  |
| ape              | 5.7.1   |
| car              | 3.1.2   |
| cowplot          | 1.1.3   |
| GET              | 1.0.4   |
| ggplot2          | 3.5.2   |
| ggpmisc          | 0.5.5   |
| ggpubr           | 0.6.0   |
| ggtext           | 0.1.2   |
| ggtree           | 3.10.1  |
| graphics         | 4.3.1   |
| ks               | 1.14.2  |
| landscapemetrics | 2.2.1   |
| mobster          | 1.0.0   |
| patchwork        | 1.2.0   |
| phangorn         | 2.11.1  |
| purrr            | 1.0.2   |
| raster           | 3.6.26  |
| spatstat         | 3.3.0   |
| spdep            | 1.3.13  |
| terra            | 1.7.78  |
| this.path        | 2.4.0   |
| tidyverse        | 2.0.0   |
| glmmTMB          | 1.1.11  |
| IOBR             | 2.2.2   |

##### Python dependencies

The simulation code was developed and tested with Python 3.x.

External Python packages required:

| Package | Version |
| :------ | :------ |
| numpy   | 1.24.3  |

The simulation scripts also use Python standard library modules:
sys, math, random, heapq, subprocess, and collections.

#### *Installation*

No compilation or installation of the source code is required.

The R analysis scripts and Python simulation scripts can be executed directly after installing the required dependencies listed above using standard package installation procedures.

The installation time depends on the existing software environment, operating system, and internet connection.

#### *Runtime*

Typical runtime for a Python simulation demo execution is less than 1 minute on a standard desktop computer. 

#### *Data source*

The dataset comprises 165 patients, including private sequencing data generated for this study (five HCC patients, DT42–DT51, and eight LUAD patients). The remaining 152 patients were collected from previously published cohorts, as listed in /data/fig1_patientID_clean.csv. For all 165 patients, VAF matrices are stored in /data/vaf_public_liver_lung_final, maximum parsimony trees are stored in /data/rds_public_liver_lung_final, and two-dimensional coordinate files are stored in /data/fig2_public_liver_lung_location_absolute_nsr_1169_final.csv.

#### *License*

The source code in this repository is licensed under the GNU General Public License v3.0. See the `LICENSE.txt` file.

Original data, metadata, documentation, and other non-software materials produced for this study are licensed under the Creative Commons Attribution 4.0 International License. See `LICENSE-DATA.txt`.

Data obtained or derived from previously published cohorts remain subject to the licensing and reuse conditions of their original sources.
