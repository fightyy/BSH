

# Deme-based Spatial Tumour Growth Simulator

This repository contains a Python script to simulate 2D/3D tumour growth and multi-region sequencing data using a deme-based agent-based model with flexible spatial constraints.

The model tracks demes (local cell population) on a 2D or 3D lattice, simulating birth, death, quiescence, mutation, deme fission, deme migration and spatial pushing. The script also generates virtual multi-region punch biopsies and corresponding VAF (variant allele frequency) matrices for downstream analysis.

## Model Overview

- **Deme-based tumour growth**  
  Each lattice site hosts a deme containing neutral and advantaged cells. Demes grow via a birth–death process until they reach a size threshold and then undergo fission.

- **Spatial constraints & pushing**  
  When no empty neighbouring site is available, deme fission can still occur via a pushing mechanism that shifts chains of demes along a direction to free space, in line with spatially constrained tumour growth models.

Requirements
---

Python version: Python **3.x**

Python packages: numpy, sys, math, random, heapq, subprocess, collections

## Usage

### 1. Set simulation dimensionality

Inside the script, the global variable controls whether you run a 2D or 3D simulation:

```bash
SIMULATION_DIM = 2  # 2 for 2D simulations, 3 for 3D simulations
```

Set this before running the script.

### 2. Command-line interface

The script is designed to be run from the command line and expects **16 arguments**:

```bash
python3 3DTumorSimulPush.py \
  <sim_dim> <deme_size> <mut_rate> <adv_rate> <s_coef> <repl> <path> <rd> \
  <birth_rate> <death_rate> <push_prop> <mig_rate> <punch_diameter> \
  <punch_density> <spacing> <title> <snapname>
```

#### Arguments

1. **`sim_dim`** (int)
    2 for 2D simulations, 3 for 3D simulations.

2. **`deme_size`** (int)
    Number of cells per deme (fission threshold is based on this).

3. **`mut_rate`** (float)
    Neutral mutation rate per cell division over the whole exonic region.

4. **`adv_rate`** (float)
    Advantageous (driver) mutation mutation rate.

5. **`s_coef`** (float)
    Selection coefficient for advantageous mutations.

6. **`repl`** (int)
    Simulation replicate index.

7. **`path`** (str)
   Output directory.

8. **`rd`** (int)

   Half-side length / radius of the lattice domain:

   - 2D: lattice size ≈ `(2*rd + 1) × (2*rd + 1)`
   - 3D: lattice size ≈ `(2*rd + 1)³`

9. **`birth_rate`** (float)
    Birth probability for cells in a deme.

10. **`death_rate`** (float)
    Death probability (must satisfy `birth_rate + death_rate ≤ 1`).
    The quiescent probability is computed as `quies_rate = 1 - birth_rate - death_rate`.

11. **`push_prop`** (float)
       Proportion of tumor diameter allowed as maximum push distance.

12. **`mig_rate`** (float)
       Migration probability per deme per event step (`MIGRATION_PROB`).

13. **`punch_diameter`** (int)
       Diameter (in lattice units) of each virtual punch biopsy used for multi-region sampling.

14. **`punch_density`** (float)
       Minimum fraction of tumour demes that must be present inside a punch for that biopsy to be considered valid.

15. **`spacing`** (int)
       Spacing (in lattice units) between the centres of adjacent punch biopsies when tiling the tumour slice.

16. **`title`** (str)
       Prefix for output file names (used for VAF files, location files, etc.).

17. **`snapname`** (str)
       Prefix for spatial snapshot files written during tumour growth.

### 3. Example command

```bash
python3 3DTumorSimulPush.py \
	2 \
  100 \
  0.3 \
  1e-5 \
  0.1 \
  1 \
  ./csv \
  50 \
  0.4 \
  0.1 \
  0.3 \
  0.05 \
  20 \
  0.95 \
  16 \
  example_sim \
  example_snap
```

This will:

- Apply deme-migration and pushing according to `mig_rate` and `push_prop`.
- Generate multi-region punch biopsies with `punch_diameter` and `spacing`.
- Write outputs to `./csv`.

------

## Outputs

All results are written to the output directory . Typical outputs include:

- **Spatial snapshots**
   Files like:
  - `<snapname>_prop1_32.txt`, …, `<snapname>_prop6_8.txt`
  - Final snapshot: `<snapname>_prop8_8.txt`
     These contain deme coordinates, deme category(if all cells in deme are neutral cells, this category will be 1 else 2 ), and mutation counts of demes.
- **Sampling locations**
  - `<title>_deme_location.txt`
     Tab-delimited file listing `(x, y, z, sample_id)` for each punch biopsy centre.
- **VAF matrix**
  - `<title>_deme_vaf.txt`
     Mutations (rows) × samples (columns), with alternating depth and VAF columns for each sample. Includes both “public” clonal mutations and private mutations.
- **Driver mutation information**
  - `<title>_adv_mutation.txt`
     Contains the ID of the advantageous mutation (`adv_id`), or `"none"` if no driver clone arose.

