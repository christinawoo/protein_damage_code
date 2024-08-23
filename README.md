# Introduction

This repository contains data on surface exposure (`relSESA`), distance from side chain amide N to backbone amide carbonyl (`distance`), and backbone torsion angles (`ramachandran_psi`, `ramachandran+phi`) for Asn and Gln residues from various datasets, as well as the code used to calculate these data. Relevant scripts are located in `/scripts` if used across multiple datasets and are otherwise co-located with the corresponding data. Scripts were either ran in UCSF ChimeraX or in a Python notebook/instance.

The repository also contains Alphafold predicted structures for truncated sequences arising from intramolecular cleavage at Asn, as well as graphics used to visualize these structures.

# Folders

## `scripts`

Global scripts used across different datasets.

## `cN_cleavage`

Analysis of Asn/Gln residues forming the C-terminal cyclic imide via intramolecular cleavage mechanism (semi-N tryptic peptides, –18.015 Da for N or Q at the peptide C-terminus). A table of cN PTM occurrences obtained from proteomics data was fed through the following scripts (in order):

1. `scripts/fetch_pdb_data.ipynb` (Python)
2. `scripts/relSESA_distance_calculation.py` (ChimeraX)
3. `scripts/ramachandran.py` (ChimeraX)
4. `scripts/json_to_excel.ipynb` (Python)
   
The resulting dataset is `cN_cleavage_data.xlsx`.

Error files generated in the calculation are located in the `errors` folder and are documented in `scripts/fetch_pdb_data.ipynb`, `scripts/relSESA_distance_calculation.py`, and `scripts/ramachandran.py`.

The AlphaFold predictions were then fed through `scripts/plddt_calc.py` (Python) to calculate the pLDDT of analyzed Asn/Gln residues and filtered for pLDDT≥70.

## `deamidation`

Analysis of Asn/Gln residues that are deamidated (fully tryptic peptides, +0.984 Da for N or Q). A table of deamidation PTM occurrences obtained from spectroscopic data was fed through the same sequence of scripts as in the `cN_cleavage` dataset.

The resulting dataset is `deamidation_data.xlsx`.

Error files generated in the calculation are located in the `errors` folder and are documented in `scripts/fetch_pdb_data.ipynb`, `scripts/relSESA_distance_calculation.py`, and `scripts/ramachandran.py`.

## `random_10000_Asn`

Analysis of Asn residues randomly selected from the [Swissprot human proteome](https://pubmed.ncbi.nlm.nih.gov/8594581/). Two sets of 5000 Asn were extracted from the human proteome using `random_10000_Asn/extract_random_asn.ipynb` (Python) and then fed through the same sequence of scripts as in the `cN_cleavage` dataset.

The resulting datasets are `random_5000_first/random_5000_first_data.xlsx` and `random_5000_second/random_5000_second_data.xlsx`, which are collated in `random_10000_ASN_combined.xlsx`.

Error files generated in the calculation are located in the `errors` folder and are documented in `scripts/fetch_pdb_data.ipynb`, `scripts/relSESA_distance_calculation.py`, and `scripts/ramachandran.py`.

## `recalculate`

Selected recalculations of relevant parameters for Asn and Gln residues from the above datasets, but using specific PDBs not initially analyzed to unify the different subunits of the same protein or substitute with more extensively cited structures. The following scripts were utilized (in order):

1. `scripts/preprocess_excel_to_json.py` (Python)
2. `scripts/analyze_specific_pdb.py` (ChimeraX)
3. `scripts/json_to_excel.ipynb` (Python)

The resulting datasets are `recalculate_first.xlsx` and `recalculate_second.xlsx`.

## `intein`

Analysis of Asn residue in Mxe GyrA intein PDBs. The following scripts were utilized (in order):

1. `intein/preprocess_analyze_intein_pdbs.py` (Python)
2. `intein/analyze_intein_pdbs.py` (ChimeraX)

The resulting dataset is `intein_cN_data.xlsx`.

## `cleaved_N_Alphafold`

Analysis of terminal Asn residues in Alphafold predicted structures of proteins truncated at the most frequent cN cleavage sites.

The dataset is `Alphafold_cleaved_N_data_pLDDT.xlsx`, which was generated using the following scripts (in order):

1. `scripts/preprocess_excel_to_json.py` (Python)
2. `cleaved_N_Alphafold/analyze_cleaved_N_Alphafold.py` (ChimeraX)
3. `scripts/json_to_excel.ipynb` (Python)

The Alphafold predictions and overlay graphics for the above dataset are in the `Alphafold_cleaved_N_predictions` folder.

Additional graphics and code used to generate graphics files are in the `graphics` folder.

### `GSS-N470 (2hgs)`

Comparison of the protein GSS (pdb `2hgs`) with the Alphafold prediction for the truncated variant (amino acids 1–470). 

Analysis of the full-length protein (pdb `2hgs`) is in `GSS_original_data.xlsx` and was obtained in a similar manner to the `recalculate` dataset.

The Alphafold prediction of the truncated variant is in `GSS_1-470_Alphafold`. The analysis of this structure is in `GSS_1-470_data.xlsx` and was obtained in a similar manner to the `cleaved_N_Alphafold` dataset.

