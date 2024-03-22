import json
import csv
import numpy as np
from chimerax.core.commands import run

def isNaN(num):
    return num != num

def runCmd(cmd):
    return run(session, cmd)

# Working directory is desktop by default: https://www.cgl.ucsf.edu/chimerax/docs/user/commands/pwd.html
data = {}

# Load preprocessed JSON
with open("/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/cleaved_n_docking/cleaved_n_docking_output_all_030124.json") as json_file:
    data = json.load(json_file)

def get_ramachandran(uniprot, is_experimental, pdb_used, chain, aa_position, original_aa):

    if isNaN(aa_position) or not aa_position:
        raise FloatingPointError('aa_position not found')

    if is_experimental:
        # Load pdb
        runCmd("open " + pdb_used)
    else:
        runCmd("alphafold fetch " + uniprot)
    
    # Renumber according to uniprot sequence
    runCmd("setattr model res_numbering uniprot")

    aa = int(aa_position)

    chain_string = "" if isNaN(chain) or not chain else  "/" + chain

    residues = runCmd('select ' + chain_string + ":" + str(aa)).residues

    ramachandran_phi = residues[0].phi
    ramachandran_psi = residues[0].psi

    switched_to_alphafold = False

    chosen_position = None

    if ramachandran_phi == None or ramachandran_psi == None:
        # Switch to alphafold
        runCmd("close all")
        runCmd("alphafold fetch " + uniprot)
        runCmd("setattr model res_numbering uniprot")

        # Select the appropriate residue
        chosen_position = str(int(original_aa))
        residue_names = runCmd("sel " + ":" +
                            chosen_position).residues.names

        # If AA position not in protein, Chimera will return an empty array
        if len(residue_names) == 0:
            raise FloatingPointError('Original AA position not found in protein after switching to alphafold: ' +
                                    uniprot + '@' + str(int(original_aa)))
        chosen_residue = residue_names[0]

        # Try the preceding residue if the given residue is not an N/Q. Accounts for methionine trimming
        if chosen_residue != 'ASN' and chosen_residue != 'GLN':
            former_position = str(int(aa_position) - 1)
            former_residue = runCmd(
                "sel " + ":" + former_position).residues.names[0]
            if former_residue == 'ASN' or former_residue == 'GLN':
                chosen_residue = former_residue
                chosen_position = former_position
            else:
                # Although not really a FPE, this error class is not used in python so we can uniquely catch this one
                raise FloatingPointError('N/Q not found in original AA or the preceding one: ' +
                                        uniprot + '@' + str(aa_position))

        residues = runCmd('select ' + ":" + chosen_position).residues

        ramachandran_phi = residues[0].phi
        ramachandran_psi = residues[0].psi

        switched_to_alphafold = True

    runCmd("close all")

    return {
        "phi": ramachandran_phi,
        "psi": ramachandran_psi,
        "switched_to_alphafold": switched_to_alphafold,
        "alphafold_switch_aa_used": chosen_position
    }

errors = []

length = len(data['uniprot_id'])

for i in map(str, range(length)):
    uniprot_id = data['uniprot_id'][i]
    is_experimental = data['is_experimental'][i]
    pdb_used = data['pdb_used'][i]
    pdb_chain = data['pdb_chain'][i]
    analyzed_position = data['analyzed_position'][i]
    original_position = data['aa_position'][i]
    try:
        ramachandran = get_ramachandran(uniprot_id, is_experimental, pdb_used, pdb_chain, analyzed_position, original_position)
        data['ramachandran_phi'][i] = ramachandran['phi']
        data['ramachandran_psi'][i] = ramachandran['psi']
        data['ramachandran_switched_to_alphafold'][i] = ramachandran['switched_to_alphafold']
        data['alphafold_switch_aa_used'][i] = ramachandran['alphafold_switch_aa_used']
    except Exception as e:
        runCmd("close all")
        position = "" if isNaN(analyzed_position) or not analyzed_position else str(int(analyzed_position))
        errors.append((uniprot_id, "Error in " + uniprot_id + '@' + position + ': ' + str(e)))

# Write output to json
with open('/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/cleaved_n_docking/cleaved_n_docking_ramachandran_030124.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

# Write errors to error file
with open('/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/cleaved_n_docking/cleaved_n_docking_ramachandran_errors_all_030124.csv', 'w') as csvfile:
    fwriter = csv.writer(csvfile)

    for x in errors:
        fwriter.writerow(x)