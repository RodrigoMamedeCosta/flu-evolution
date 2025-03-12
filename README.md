# **Rolling Sphere and Statistical Analyses**

## **Overview**
This repository contains **R scripts and datasets** required for reproducing the rolling sphere analysis and most statistical analyses from:

> **Costa RM, Acosta-Alvarez L, Curtis K, Zarbock K, Kelleher J, et al.**  
> *Influenza virus evolution and defective genome formation are shaped by host genotype and sex* (2025).  
> [[10.1101/2025.02.26.638946v1](https://www.biorxiv.org/content/10.1101/2025.02.26.638946v1)]

## **Contents**
- `3D structures with replaced b-factors/` – These are the modified PDBs for chimera/chimerax, used to generate the visualizations. (Fig. 4)

- `Deletion analysis/` – This contains the unfiltered read counts (ancestral and derived) and the script that looks for deletions based
   on contiguous stretches of derived reads at similar frequency. Also includes scripts for statistical analyses of deletions. (Fig. 6)

- `Rolling sphere/` – Script for the rolling sphere model used to detect positive selection, along with the required data (in the format of
   filtered ancestral and derived reads and PDB-derived distance matrices).

- `Titer Analysis/` – Contains data and scripts for statistical analysis of differences in viral titers (fig. 2)

- `Weight loss and survival/` – Contains data and scripts for statistical analysis of differences in virulence (fig. 2, 3, 7)

## **Installation & Requirements**
Ensure you have **R 4.0+** installed. Required packages are included in scripts.

## Licensed under the **MIT License**

## **Citation**
If you use these scripts in your research, please cite:

> **Costa RM, Acosta-Alvarez L, Curtis K, Zarbock K, Kelleher J, et al.**  
> *Influenza virus evolution and defective genome formation are shaped by host genotype and sex* (2025).  
> [10.1101/2025.02.26.638946v1]

## **Contact**
For questions or collaborations, reach out to:
- Rodrigo Costa
- rodrigomcosta@gmail.com
- University of Utah
