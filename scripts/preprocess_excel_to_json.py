"""
  Script for converting .xlsx files to .json files, primarily for use as input to analyze_specific_pdb.py or analyze_cleaved_N_Alphafold.py.

  
  For use with analyze_specific_pdb.py, the Excel spreadsheet should contain the following columns:
    uniprot_id: Required. Uniprot accession of protein to be analyzed
    aa_position: Required. Position of Asn/Gln of interest within protein
    pdb_id: Required. PDB code to use
    pdb_chain: Required. PDB chain to use

    analyzed_position: Blank
    analyzed_aa: Blank
    calculated_SES: Blank
    relSESA: Blank
    distance: Blank
    ramachandran_phi: Blank
    ramachandran_psi: Blank

    
  For use with analyze_cleaved_N_Alphafold.py, the Excel spreadsheet should contain the following columns:

    gene_name_and_position: Required. Concatenated Uniprot accession of protein and the position of Asn/Gln to be analyzed. E.g., for Uniprot accession P68871 and Asn position 58, gene_name_and_position should be denoted as P6887158.
        
    analyzed_aa_position: Required. Position of Asn/Gln of interest within protein
    pdb_chain: Required. Chain in PDB to use (usually A)

    analyzed_aa: Blank
    calculated_SES: Blank
    relSESA: Blank
    distance: Blank
    ramachandran_phi: Blank
    ramachandran_psi: Blank
"""

INPUT_EXCEL_FILEPATH = "" # (Absolute) path to input .xlsx file (include .xlsx extension)
OUTPUT_JSON_FILEPATH = "" # (Absolute) path to output .json file (include .json extension)

import pandas as pd
import json

df = pd.read_excel(INPUT_EXCEL_FILEPATH)

df_dict = df.to_dict()

with open(OUTPUT_JSON_FILEPATH, 'w', encoding='utf-8') as f:
    json.dump(df_dict, f, ensure_ascii=False, indent=4)
