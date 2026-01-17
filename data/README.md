### Folder structure

- `vaf_public_liver_lung/`  
  Variant allele frequency (VAF) matrices for 165 patients.  
  
  - Filename: `tumortype_patientID`  
  - `mut_id`: chromosome, position, reference base, and alternate base joined by `_` (e.g. `chr1_12345_A_T`)  
  - Remaining columns: sample names; values = mutation frequency (VAF)
  
- `rds_public_liver_lung/`  
  Maximum parsimony trees for the same 165 patients.  
  
  - Filename: `tumortype_patientID`  
  - Each file is an RDS object storing the patient-level tree.
  
- `fig2_monmonier_out_two_path/`  
  Boundary coordinates for 23 patients with ≥10 samples.  
  
  - `bound_x`: x-coordinate of boundary  
  - `bound_y`: y-coordinate of boundary  
  - `run_id`: boundary index  
  - `dir`: boundary direction  
  - `Step`: step index of edges composing this boundary  
  - `patient_cancer`: `tumortype_patientID`  
  - `tumor_type`: tumor type  
  - `threshold`: when an edge is classified as boundary, the lowest quantile of the genetic distance of this edges  
  - `sample_num`: number of samples for this patient  
  - ###### `csv_example/`
  
    Example tumor data used for plotting. Filenames follow the conventions and column definitions of `fig3_hu_result_hu_time_result.csv`.
  
    - Files ending with `location`  
      - Example tables for sample locations  
      - `sample`: sample name  
      - `x`, `y`, `z`: center coordinates of the bulk sample  
  
    - Files ending with `vaf`  
      - Example VAF tables  
      - `mut_id`: mutation ID  
      - `depth1`: sequencing depth of each mutation in sample 1  
      - `S1`: mutation frequency (VAF) of each mutation in sample 1  
  
    - Files ending with `boundary`  
      - Boundary strength files  
      - Column definitions: see `fig3_hu_time_boundary_silhouette_freq.csv`  
  
    - Files ending with `result`  
      - Simulated ENN values  
      - Column definitions: see `fig3_hu_result_hu_time_result.csv`  
  
    - Files ending with `tree`  
      - Maximum parsimony trees  
  
    - Files ending with `prop8_8`  
      - Coordinates of the central tumor slice  
      - `x, y`: coordinates  
      - `category`: if `category == 2`, the deme contains selected (advantageous) cells  
      - `mutation_count`: number of mutations contained in this deme  

---

### Other CSV files

### Figure 1

#### `fig1_patientID_clean.csv`

This file contains metadata for all 165 patients. Column definitions are as follows:

- `Abbreviation`: Abbreviation of the tumor type.
- `PatientID`: Unique identifier assigned to each patient.
- `Tumor type`: Tumor type.
- `Patient`: Original patient identifier in the source study.
- `Sample number`: Number of samples collected for this patient.
- `Paper / Data Collection Name`: Name of the source publication or data collection.
- `DOI`: DOI of the source publication.
- `Year`: Publication year of the source study.

#### `fig1_result_slope_all_relative_summary.csv`  

Summary of linear regressions between genetic distance and relative physical distance for all patients. Detailed pairwise data: `fig1_result_slope_all_relative.csv`.  
- `patient`: `tumortype_patientID`  
- `slope`: slope of the patient-level linear regression  
- `slope_p_value`: p-value for the slope  
- `intercept`: intercept of the linear regression  
- `intercept_p_value`: p-value for the intercept  
- `r_squared`: R2 of the regression  
- `tumor_type`: tumor type  
- `sample_number`: number of samples for this patient  

#### `fig1_result_slope_all_relative.csv`  
Relative physical distances and genetic distances for all patient–sample pairs.  
- `Sample1`: first sample in the pairwise comparison  
- `Sample2`: second sample in the pairwise comparison  
- `ITH_values`: (private mutations in Sample1 + private mutations in Sample2) / (all mutations in Sample1 and Sample2), range 0–1  
- `Physical_Dist`: physical distance between Sample1 and Sample2  
- `Tumor_type`: tumor type  
- `Patient`: `tumortype_patientID`  
- `max_dist`: maximum distance between any two samples within this patient  
- `Relative_physical_distance`: `Physical_Dist / max_dist`  

### Figure 2

#### `fig2_public_liver_lung_loaction_absolute_nsr.csv`  
2D coordinates of all samples for 165 patients.  
- `Patient`: patient ID  
- `Sample`: sample name  
- `Tumor_type`: tumor type  
- `X`: x-coordinate  
- `Y`: y-coordinate  

#### `fig2_result_boundary_real.csv`  
Boundary strength between genetic blocks for all patients (real data). Only boundaries between spatially neighboring genetic classes are considered.  
- `tumor_type`: tumor type  
- `Patient`: `tumortype_patientID`  
- `cat1`: first genetic class in the pairwise comparison  
- `cat2`: second genetic class in the pairwise comparison  
- `n_cat1`: number of samples in `cat1`  
- `n_cat2`: number of samples in `cat2`  
- `silhouette_mean`: silhouette value between `cat1` and `cat2`, representing boundary strength  

#### `fig2_result_enn_treestat_real_data.csv`  
Blockiness values ENNmn for all patients (real data).  
- `tumor_type`: tumor type  
- `Patient`: `tumortype_patientID`  
- `value_enn`: ENNmn value for each patient; ENNmn = 0 indicates strictly blocky structure, larger ENNmn indicates weaker blockiness  

---

### Figure 3 

#### `fig3_hu_result_hu_time_result.csv`  
Simulated ENNmn values under different selection coefficients and push power settings.  
- `s_coef`: selection coefficient; birth rate = `(1 + s_coef) * 0.4`  
- `adv_rate`: advantageous mutation rate  
- `push_prop`: push power, proportion of tumor diameter allowed as maximum push distance  
- `mut_rate`: mutation rate (mutations per cell division)  
- `sample_diameter`: diameter of each virtual punch biopsy (in demes) used for multi-region sampling  
- `vaf_cutoff`: VAF threshold to call a mutation  
- `seed`: random seed for the simulation  
- `diameter`: tumor diameter (in demes)  
- `value_enn`: ENNmn value for each simulation; ENNmn = 0 indicates strictly blocky structure, larger ENNmn indicates weaker blockiness  

#### `fig3_hu_time_boundary_silhouette_freq.csv`  
Boundary strength of simulated tumors under different selection coefficients.  
- `s_coef`: selection coefficient; birth rate = `(1 + s_coef) * 0.4`  
- `adv_rate`: advantageous mutation rate  
- `seed`: random seed for the simulation  
- `diameter`: tumor diameter (in demes)  
- `cat1`: first genetic class in the pairwise comparison  
- `cat2`: second genetic class in the pairwise comparison  
- `n_cat1`: number of samples in `cat1`  
- `n_cat2`: number of samples in `cat2`  
- `silhouette_mean`: mean silhouette between `cat1` and `cat2`, representing boundary strength  
- `whole_adv_freq`: mutation frequency of advantageous mutation across the tumor  

#### `fig3_hu_time_boundary_silhouette_push.csv`  
Boundary strength of simulated tumors under different push power values.  
- Column definitions: same as `fig3_hu_time_boundary_silhouette_freq.csv`  
- `push_prop`: push power, proportion of tumor diameter allowed as maximum push distance  

#### `fig3_hu_time_boundary_silhouette.csv`  
Boundary strength of neutral simulations (no selection).  
- Column definitions: same as `fig3_hu_time_boundary_silhouette_freq.csv`  

---

### Figure 4

#### `fig4_go_enrich.csv`  
GO enrichment results for differentially expressed genes between HCC tumors with high ENNmn and low ENNmn.  
- `ID`: GO ID  
- `Description`: short description of the GO biological process / pathway / function  
- `GeneRatio`: proportion of genes in the input gene list associated with this term  
- `BgRatio`: proportion of background genes associated with this term  
- `pvalue`: raw p-value for enrichment  
- `p.adjust`: p-value after multiple-testing correction (Benjamini–Hochberg)  
- `qvalue`: estimated false discovery rate (FDR)  
- `geneID`: input genes associated with this term, separated by `/`  
- `Count`: number of input genes falling into this term  
- `Direction`: direction of regulation for the gene set (e.g. `Up` or `Down`)  
- `Score`: enrichment score,  `−log10(qvalue)`  

#### `fig4_gsva_enn.csv`  
GSVA scores for GO:0098742 in HCC patients and corresponding ENNmn values.  
- `sample`: sample name, formatted as `patientID_sampleID`  
- `GO_0098742_score`: GSVA score for GO:0098742; >0 indicates up-regulation, <0 indicates down-regulation  
- `value_enn`: ENNmn value for the corresponding patient  

#### `fig4_hu_time_boundary_silhouette_migration.csv`  
Boundary strength of simulated tumors under different migration rates.  
- Column definitions: same as `fig3_hu_time_boundary_silhouette_freq.csv`  
- `migration rate`: migration probability per deme per event step  

#### `fig4_real_data_center_edge_ITH.csv`  
ITH values of center and edge regions for all patients (real data).  
- `location`: center or edge; a sample is labeled edge if its distance to the nearest tumor boundary is <25% of the maximum tumor diameter, otherwise center
- `ITH_mean`: mean ITH in the center or edge region. ITH is defined as the average genetic distance between all pairs of samples within a circular area (radius = median pairwise distance of samples in that area).  

#### `fig4_simulation_center_edge_ITH.csv`  
ITH values of center and edge regions in boundary-growth simulations.  
- `ITH_center`: mean ITH in the center region  
- `ITH_edge`: mean ITH in the edge region  

#### `fig4_result_center_edge_real_simulation.csv`  
Mean mutation counts of center and edge samples in boundary-growth simulations and real data.  
- `edge`: mean number of mutations in edge samples (center/edge defined as above using distance to nearest tumor boundary)  
- `center`: mean number of mutations in center samples  

#### `fig4_result_enn_treestat_migrate.csv`  
ENN values under different migration rates.  
- Column definitions: see `fig3_hu_result_hu_time_result.csv`  
- `migration rate`: migration probability per deme per event step  

---

### Figure 5

#### `fig5_hu_time_boundary_mob_selection_single_birth0.4_0.1-0.9.csv`  
MOBSTER output per sample in simulations with selection.  
- General parameter columns: see `fig3_hu_result_hu_time_result.csv`  
- `mobster_result_vec`: MOBSTER result for each sample  
- `sample_type_vec`: genetic class of each sample  
- `x`: x-coordinate of the sample  
- `y`: y-coordinate of the sample  
- `sample`: sample name  
- `mobster_test_type`: MOBSTER test type (`single` for single samples, `window` for sliding windows)  
- `whole_tumor_freq`: mutation frequency of advantageous mutation in the whole tumor  
- `Tail`: number of mutations in the MOBSTER neutral tail  
- `C1`, `C2`, `C3`: number of mutations in MOBSTER clusters 1, 2, and 3  
- `Boundary`: whether the sample lies on a boundary between genetic classes  

#### `fig5_hu_time_boundary_mob_selection_window_birth0.4_0.1-0.9.csv`  
MOBSTER output per window in simulations with selection.  
- Column definitions: see `fig5_hu_time_boundary_mob_selection_single_birth0.4_0.1-0.9.csv`  

#### `fig5_mob_boundary_real_data.csv`  
MOBSTER results for each sample in real data.  
- `Punch`: sample name, typically `patientID_sampleID`  
- `Boundary`: whether the sample lies on a genetic boundary  
- `mobster_test_type`: MOBSTER test type (`single` or `window`)  
- `mobster_result`: MOBSTER result for each sample  
- `genetic_type`: genetic class of the sample  
- `input_type`: type of input VAF, including raw VAF, mpileup VAF, CCF maximum a posteriori estimate, and CCF posterior mean  

#### `fig5_mob_pileup_result_ccf_vaf_window_single.csv`  
Sliding-window MOBSTER results for real data.  
- Column definitions: see `fig5_mob_boundary_real_data.csv`  

#### `fig5_mob_single_window_neutral.csv`  
MOBSTER output for each sample and sliding window in neutral simulations.  
- Column definitions: see `fig5_hu_time_boundary_mob_selection_single_birth0.4_0.1-0.9.csv`  

#### `fig5_result_center_edge_all0.25.csv`  
Center vs edge classification for HCC and LUAD samples.  
- `location`: center or edge, defined by distance to the nearest tumor boundary; samples with distance <25% of the maximum tumor diameter are labeled edge, otherwise center  

#### `fig5_result_dnds_singles_sample.csv`  
dN/dS values per sample for HCC and LUAD.  
- `name`: dN/dS gene category, including missense (`wmis`), nonsense (`wnon`), essential splice site (`wspl`), indels (`wind`), all non-synonymous substitutions (`wall`), and truncating substitutions (`wtru`)  
- `mle`: maximum-likelihood estimate  
- `cilow`: lower bound of the confidence interval  
- `cihigh`: upper bound of the confidence interval  
- `SampleID`: sample name  

#### `fig5_result_sim_center_edge_neutral_0.25.csv`

Center vs edge classification for neutral simulated tumor.

- `location`: center or edge, defined by distance to the nearest tumor boundary; samples with distance <25% of the maximum tumor diameter are labeled edge, otherwise center  