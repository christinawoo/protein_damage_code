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
    # Delete solvent
    runCmd("delete solvent")
    if isinstance(chain, str):
      # Remove other chains
      runCmd("sel /" + chain)
      runCmd("select ~sel")
      runCmd("delete sel")
  else:
    # Fetch alphafold structure
    runCmd("alphafold fetch " + uniprot)

  # Set chain identifier depending on if chain exists
  chain_id = "/" + chain if isinstance(chain, str) else ""

  # Select the appropriate residue
  chosen_position = str(aa_position)
  residue_names = runCmd("sel " + chain_id + ":" + chosen_position).residues.names
  if len(residue_names) == 0:
    raise FloatingPointError('AA position not found in protein: ' + uniprot + ', ' + chain_id + '@' + str(aa_position))
  chosen_residue = residue_names[0]

  # Try the preceding residue if the given residue is not an N/Q. Accounts for methionine trimming
  if chosen_residue != 'ASN' and chosen_residue != 'GLN':
    former_position = str(int(aa_position) - 1)
    former_residue = runCmd("sel " + chain_id + ":" + former_position).residues.names[0]
    if former_residue == 'ASN' or former_residue == 'GLN':
      chosen_residue = former_residue
      chosen_position = former_position
    else:
      # Although not really a FPE, this error class is not used in python so we can uniquely catch this one
      raise FloatingPointError('N/Q not found in this AA or the preceding one: ' + uniprot + ', ' + chain_id + '@' + str(aa_position))
  
  # Calculate distance
  side_chain_N = chain_id + ':' + chosen_position + '@ND2' if chosen_residue == 'ASN' else chain_id + ':' + chosen_position + '@NE2'
  carbonyl_C = chain_id + ':' + chosen_position + '@C'

  print(side_chain_N)
  print(carbonyl_C)
  command = 'distance ' + side_chain_N + ' ' + carbonyl_C
  print(command)
  distance = runCmd(command)

  # Calculate SESA
  runCmd("sel " + chain_id + ":" + chosen_position)
  runCmd("surface sel")
  ses = runCmd("measure area sel includeMasked false")

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

errors = []

length = len(data['gene_name_and_position'])
for i in map(str, range(length)):


  try:
    uniprot_id = data['uniprot_id'][i]
    is_experimental = data['is_experimental'][i]
    pdb = data['pdb_used'][i]
    chain = data['pdb_chain'][i]
    aa_position = data['aa_position'][i]

    analysis = analyze(uniprot_id, is_experimental, pdb, chain, aa_position)

    # Any errors will be caught above 
    data['analyzed_position'][i] = analysis['analyzed_position']
    data['analyzed_aa'][i] = analysis['aa']
    data['calculated_SES'][i] = analysis['ses']
    data['relSESA'][i] = analysis['relSESA']
    data['distance'][i] = analysis['distance']

  except Exception as e:
    # Important to close the running file if there is an error!
    runCmd('close all')
    chain_id = "/" + chain if isinstance(chain, str) else ""
    errors.append((uniprot_id, "Error in " + uniprot_id + ', ' + chain_id + '@' + str(aa_position) + ': ' + str(e)))

# Write output to json
with open('./chimera/output_all.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

# Write errors to error file
with open('./chimera/errors_all.csv', 'w') as csvfile:
    fwriter = csv.writer(csvfile)

    for x in errors:
        fwriter.writerow(x)

