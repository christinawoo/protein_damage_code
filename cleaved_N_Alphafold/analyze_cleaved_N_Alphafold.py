"""
    Script for calculating relative solvent excluded surface area (relSESA), distance from side chain N to backbone carbonyl C, and backbone torsion angles for Asn and Gln residues in Alphafold predicted models of proteins truncated at the cN cleavage site.

    The input PDBs should be generated via Alphafold and must be organized in a folder with the following structure:

    - Root folder (e.g. "Alphafold_cleaved_N_predictions")
        - {gene_name_and_position}
            - {gene_name_and_position}_selected_prediction.pdb
    
    where gene_name_and_position is the concatenated Uniprot accession of protein and the position of Asn/Gln to be analyzed. E.g., for Uniprot accession P68871 and Asn position 58, gene_name_and_position should be denoted as P6887158, and the PDB should be located at Alphafold_cleaved_N_predictions/P6887158/P6887158_selected_prediction.pdb

    The input to this script should be a JSON file obtained from preprocess_excel_to_json.py, which should accept an Excel file containing the following columns:
        
        gene_name_and_position: Required. See description above.
        analyzed_aa_position: Required. Position of Asn/Gln of interest within protein
        pdb_chain: Required. Chain in PDB to use (usually A)

        analyzed_aa: Blank
        calculated_SES: Blank
        relSESA: Blank
        distance: Blank
        ramachandran_phi: Blank
        ramachandran_psi: Blank
"""

# Absolute path to root folder (e.g. Alphafold_cleaved_N_predictions) containing input PDBs (see description above).
PDB_ROOT_FOLDER_PATH = ""

# Path to input file (should include the desired filename ending in .json)
INPUT_JSON_FILEPATH = ""

# Path to output JSON file (should include the desired filename ending in .json)
OUTPUT_JSON_FILEPATH = ""

import json
import csv
import numpy as np
from chimerax.core.commands import run

def runCmd(cmd):
    return run(session, cmd)

def analyze(gene_name_and_position, chain, aa_position):
    print(gene_name_and_position)
    print(aa_position)
    # Load structure
    # Load pdb
    runCmd(f'open "{PDB_ROOT_FOLDER_PATH}/{gene_name_and_position}/{gene_name_and_position}_selected_prediction.pdb"')

    # ! Removed this for cleaved n docking, as no cleanup is needed.
    # # Remove other chains
    # runCmd("sel /" + chain)
    # runCmd("select ~sel")
    # runCmd("delete sel")

    # # Delete solvent
    # runCmd("delete solvent")

    # # Delete ligand
    # runCmd("delete ligand")

    # ! Removed this for cleaved n docking to avoid renumbering errors, as exact AA positions to be analyzed in the PDB are known.
    # Renumber according to uniprot sequence
    # runCmd("setattr model res_numbering uniprot")

    # Set chain identifier
    chain_id = "/" + chain

    # Select the appropriate residue
    chosen_position = str(aa_position)
    residue_names = runCmd("sel " + chain_id + ":" +
                           chosen_position).residues.names

    # If AA position not in protein, Chimera will return an empty array
    if len(residue_names) == 0:
        raise FloatingPointError('AA position not found in protein: ' +
                                 gene_name_and_position + ', ' + chain_id + '@' + str(aa_position))
    chosen_residue = residue_names[0]

    # ! Commented out the below check for preceding residues for cleaved N docking, because the exact residue position is already known
    # Try the preceding residue if the given residue is not an N/Q. Accounts for methionine trimming
    # if chosen_residue != 'ASN' and chosen_residue != 'GLN':
    #     former_position = str(int(aa_position) - 1)
    #     former_residue = runCmd(
    #         "sel " + chain_id + ":" + former_position).residues.names[0]
    #     if former_residue == 'ASN' or former_residue == 'GLN':
    #         chosen_residue = former_residue
    #         chosen_position = former_position
    #     else:
    #         # Although not really a FPE, this error class is not used in python so we can uniquely catch this one
    #         raise FloatingPointError('N/Q not found in this AA or the preceding one: ' +
    #                                  uniprot + ', ' + chain_id + '@' + str(aa_position))

    # Calculate distance
    side_chain_N = chain_id + ':' + chosen_position + \
        '@ND2' if chosen_residue == 'ASN' else chain_id + ':' + chosen_position + '@NE2'
    carbonyl_C = chain_id + ':' + chosen_position + '@C'
    command = 'distance ' + side_chain_N + ' ' + carbonyl_C
    distance = runCmd(command)
    print(side_chain_N)
    print(carbonyl_C)
    print(command)

    # Calculate SESA
    runCmd("sel " + chain_id + ":" + chosen_position)
    runCmd("surface sel")
    ses = runCmd("measure area sel includeMasked false")

    # Assumes the residue iseither ASN or GLN, values from areaSESgxg.txt
    relSESA = ses / 90.541 if chosen_residue == 'ASN' else ses / 106.534

    # Calculate Ramachandran
    # ! Need to use the analyzed position, not the original position, to account for off by one errors!
    aa = int(chosen_position)

    chain_string = "/" + chain

    residues = runCmd('select ' + chain_string + ":" + str(aa)).residues

    ramachandran_phi = residues[0].phi
    ramachandran_psi = residues[0].psi

    # Close
    runCmd('close all')

    # Return values
    output = {
        "analyzed_position": chosen_position,
        "aa": chosen_residue,
        "distance": distance,
        "ses": ses,
        "relSESA": relSESA,
        "phi": ramachandran_phi,
        "psi": ramachandran_psi,
    }
    return output


# Working directory is desktop by default: https://www.cgl.ucsf.edu/chimerax/docs/user/commands/pwd.html
data = {}

# Load preprocessed JSON
with open(INPUT_JSON_FILEPATH) as json_file:
    data = json.load(json_file)

# Number of data points
length = len(data['gene_name_and_position'])

for i in map(str, range(length)):

    # Define the uniprot id and aa position early
    gene_name_and_position = data['gene_name_and_position'][i]
    aa_position = data['analyzed_aa_position'][i]
    chain = data['pdb_chain'][i]

    analysis = analyze(gene_name_and_position, chain, aa_position)

    # Any errors will be caught above
    # ! Removed this because you the position to be analyzed is already known
    # data['analyzed_position'][i] = analysis['analyzed_position']
    data['analyzed_aa'][i] = analysis['aa']
    data['calculated_SES'][i] = analysis['ses']
    data['relSESA'][i] = analysis['relSESA']
    data['distance'][i] = analysis['distance']
    data['ramachandran_phi'][i] = analysis['phi']
    data['ramachandran_psi'][i] = analysis['psi']

# Write output to json
with open(OUTPUT_JSON_FILEPATH, 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
