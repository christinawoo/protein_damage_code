import json
import csv
import numpy as np
from chimerax.core.commands import run


def runCmd(cmd):
    return run(session, cmd)


def analyze(uniprot, is_experimental, pdb, chain, aa_position):
    print(uniprot)
    # Load structure
    if is_experimental:
        # Load pdb
        runCmd("open " + pdb)

        # If there is a chain (string) to select, remove other chains
        if isinstance(chain, str):
            # Remove other chains
            runCmd("sel /" + chain)
            runCmd("select ~sel")
            runCmd("delete sel")
    else:
        # Fetch alphafold structure
        runCmd("alphafold fetch " + uniprot)

    # Delete solvent
    runCmd("delete solvent")

    # Delete ligand
    runCmd("delete ligand")

    # Renumber according to uniprot sequence
    runCmd("setattr model res_numbering uniprot")

    # Set chain identifier depending on if chain exists
    chain_id = "/" + chain if isinstance(chain, str) else ""

    # Select the appropriate residue
    chosen_position = str(aa_position)
    residue_names = runCmd("sel " + chain_id + ":" +
                           chosen_position).residues.names

    # If AA position not in protein, Chimera will return an empty array
    if len(residue_names) == 0:
        raise FloatingPointError('AA position not found in protein: ' +
                                 uniprot + ', ' + chain_id + '@' + str(aa_position))
    chosen_residue = residue_names[0]

    # Try the preceding residue if the given residue is not an N/Q. Accounts for methionine trimming
    if chosen_residue != 'ASN' and chosen_residue != 'GLN':
        former_position = str(int(aa_position) - 1)
        former_residue = runCmd(
            "sel " + chain_id + ":" + former_position).residues.names[0]
        if former_residue == 'ASN' or former_residue == 'GLN':
            chosen_residue = former_residue
            chosen_position = former_position
        else:
            # Although not really a FPE, this error class is not used in python so we can uniquely catch this one
            raise FloatingPointError('N/Q not found in this AA or the preceding one: ' +
                                     uniprot + ', ' + chain_id + '@' + str(aa_position))

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

    # Close
    runCmd('close all')

    # Return values
    output = {
        "analyzed_position": chosen_position,
        "aa": chosen_residue,
        "distance": distance,
        "ses": ses,
        "relSESA": relSESA
    }
    return output


# Working directory is desktop by default: https://www.cgl.ucsf.edu/chimerax/docs/user/commands/pwd.html
data = {}

# Load preprocessed JSON
with open("./chimera/preprocess.json") as json_file:
    data = json.load(json_file)

# Errors occuring in individual tries of PDBs
errors = []

# Errors caused when all PDBs for a row have been exhausted, thus resorting to alphafold.
# This is different from a row not having available PDB's in the first place (i.e. alphafold by default)
all_failed = []

# Errors in Alphafold computations. This should be empty
alphafold_errors = []

# Number of data points
length = len(data['gene_name_and_position'])

# ! This should be range length in future iterations
for i in map(str, range(2000)):

    # How many PDBs there were originally for the data point
    initial_obj_length = len(data['pdb_data_obj'][i])

    # Define the uniprot id and aa position early
    uniprot_id = data['uniprot_id'][i]
    aa_position = data['aa_position'][i]

    # Process of running alphafold computation
    def try_alphafold():
        try:
            analysis = analyze(uniprot_id, False, None, None, aa_position)

            # Populate dict with the obtained analysis data
            data['analyzed_position'][i] = analysis['analyzed_position']
            data['analyzed_aa'][i] = analysis['aa']
            data['calculated_SES'][i] = analysis['ses']
            data['relSESA'][i] = analysis['relSESA']
            data['distance'][i] = analysis['distance']
            data['is_experimental'][i] = False
        except Exception as e:
            # Must close all to make sure the cycle fully terminates
            runCmd('close all')
            chain_id = ""
            alphafold_errors.append((uniprot_id, "ALPHAFOLD Error in " + uniprot_id +
                                     ', ' + chain_id + '@' + str(aa_position) + ': ' + str(e)))

    # If PDB data could not be fetched in the API process on Deepnote
    if data['pdb_data_obj'][i] == None:
        try_alphafold()

    # If PDB data could be fetched on Deepnote, but none met the criteria for resolution and range
    elif len(data['pdb_data_obj'][i]) == 0:
        try_alphafold()

    # Otherwise, attempt PDB analysis
    else:
        # List of failed pdb codes
        failed_pdbs = []

        # Repeat until all PDB objects in the object list have been exhausted
        while len(data['pdb_data_obj'][i]) > 0:

            # Current PDB data
            pdb_obj = data['pdb_data_obj'][i][0]

            # If PDB has already been tried, remove from list and finish this loop iteration.
            if pdb_obj['pdb_id'] in failed_pdbs:
                data['pdb_data_obj'][i].pop(0)
                continue

            try:
                # Set input values for analysis
                is_experimental = True
                pdb_used = pdb_obj['pdb_id']
                pdb_chain = pdb_obj['chain_id']
                resolution = pdb_obj['resolution']
                pdb_experimental_method = pdb_obj['experimental_method']

                pdb = pdb_used
                chain = pdb_chain

                analysis = analyze(uniprot_id, is_experimental,
                                   pdb, chain, aa_position)

                # Any errors will be caught above
                data['analyzed_position'][i] = analysis['analyzed_position']
                data['analyzed_aa'][i] = analysis['aa']
                data['calculated_SES'][i] = analysis['ses']
                data['relSESA'][i] = analysis['relSESA']
                data['distance'][i] = analysis['distance']

                data['is_experimental'][i] = is_experimental
                data['pdb_used'][i] = pdb_used
                data['pdb_chain'][i] = pdb_chain
                data['pdb_resolution'][i] = resolution
                data['pdb_experimental_method'][i] = pdb_experimental_method

                break
            except Exception as e:
                # Important to close the running file if there is an error!
                runCmd('close all')

                # Log error
                chain_id = "/" + chain if isinstance(chain, str) else ""
                errors.append((uniprot_id, "Error in " + pdb_obj['pdb_id'] + ", " + uniprot_id +
                              ', ' + chain_id + '@' + str(aa_position) + ': ' + str(e) + '. Skipped all other attempts on this PDB file.'))

                # Remove the failed item from the object list and add to list of failed pdb codes
                failed = data['pdb_data_obj'][i].pop(0)
                failed_pdbs.append(failed['pdb_id'])

        # This should only run if all PDBs have been exhausted
        if len(data['pdb_data_obj'][i]) == 0:
            all_failed.append((uniprot_id, "ALL FAILED Error in " + uniprot_id +
                              ', ' + chain_id + '@' + str(aa_position) + ", tried " + str(initial_obj_length) + " files"))
            try_alphafold()


# Write output to json
with open('./chimera/1-2000output_all.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

# Write errors to error file
with open('./chimera/1-2000errors_all.csv', 'w') as csvfile:
    fwriter = csv.writer(csvfile)

    for x in errors:
        fwriter.writerow(x)


with open('./chimera/1-2000errors_alphafold.csv', 'w') as csvfile:
    fwriter = csv.writer(csvfile)

    for x in alphafold_errors:
        fwriter.writerow(x)

with open('./chimera/1-2000errors_all_failed.csv', 'w') as csvfile:
    fwriter = csv.writer(csvfile)

    for x in all_failed:
        fwriter.writerow(x)
