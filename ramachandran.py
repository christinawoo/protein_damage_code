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
with open("./chimera/new_ramachandran_preprocess.json") as json_file:
    data = json.load(json_file)

def get_ramachandran(uniprot, is_experimental, pdb_used, chain, aa_position):

    if isNaN(aa_position):
        raise FloatingPointError('aa_position not found')

    if is_experimental:
        # Load pdb
        runCmd("open " + pdb_used)
    else:
        runCmd("alphafold fetch " + uniprot)
    
    # Renumber according to uniprot sequence
    runCmd("setattr model res_numbering uniprot")

    aa = int(aa_position)

    chain_string = "" if isNaN(chain) else  "/" + chain

    residues = runCmd('select ' + chain_string + ":" + str(aa)).residues

    ramachandran_phi = residues[0].phi
    ramachandran_psi = residues[0].psi

    # print(test_phi)
    # print(test_psi)

    # print('yeet')
    # print(residues)

    # ramachandran_phi = runCmd("torsion " + chain_string + ":" + str(aa) + "@c,ca,n:" + str(aa - 1) + "@c")

    # ramachandran_psi = runCmd("torsion " + chain_string + ":" + str(aa) + "@n,ca,c:" + str(aa + 1) + "@n")

    runCmd("close all")

    return {
        "phi": ramachandran_phi,
        "psi": ramachandran_psi
    }

errors = []

for i in map(str, range(2000)):
    uniprot_id = data['uniprot_id'][i]
    is_experimental = data['is_experimental'][i]
    pdb_used = data['pdb_used'][i]
    pdb_chain = data['pdb_chain'][i]
    analyzed_position = data['analyzed_position'][i]
    try:
        ramachandran = get_ramachandran(uniprot_id, is_experimental, pdb_used, pdb_chain, analyzed_position)
        data['ramachandran_phi'][i] = ramachandran['phi']
        data['ramachandran_psi'][i] = ramachandran['psi']
    except Exception as e:
        runCmd("close all")
        position = "" if isNaN(analyzed_position) else str(int(analyzed_position))
        errors.append((uniprot_id, "Error in " + uniprot_id + '@' + position + ': ' + str(e)))

# Write output to json
with open('./chimera/1-2000ramachandran_all.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

# Write errors to error file
with open('./chimera/1-2000ramachandran_errors_all.csv', 'w') as csvfile:
    fwriter = csv.writer(csvfile)

    for x in errors:
        fwriter.writerow(x)