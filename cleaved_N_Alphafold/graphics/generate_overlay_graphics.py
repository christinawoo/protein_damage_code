"""
    Script used to generate consistent overlays of cleaved-N Alphafold predictions with the original full-length proteins.

    The input to this script should be a JSON file obtained from preprocess_excel_to_json.py, which should accept an Excel file containing the following columns:
        gene_name_and_position
        original_position
        cleaved_position 
        chain 
        pdb 
        label
        pdb_relSESA 
        pdb_distance
        truncated_relSESA
        truncated_distance

    The Excel file used in this study is graphics_generation_data.xlsx.
"""

# Path to input file (should include the desired filename ending in .json)
INPUT_JSON_FILEPATH = ""

# Set row in input JSON to use to generate graphic. Feel free to change this!
row = str(1)

from chimerax.core.commands import run
import json

# ! Initially run 'windowsize 1000 600' to standardize window size
# ! Before exporting, adjust zoom with 'zoom pixelSize 0.02'


# Load preprocessed JSON
with open(INPUT_JSON_FILEPATH) as json_file:
    data = json.load(json_file)

gene_name_and_position = data["gene_name_and_position"][row]
original_position = data["pdb_position"][row]
cleaved_position = data["cleaved_position"][row]
chain = data["pdb_chain"][row]
pdb = data["pdb"][row]
label = data["label"][row]

pdb_relSESA = data["original_relSESA"][row]
pdb_distance = data["original_distance"][row]

truncated_relSESA = data["cleaved_relSESA"][row]
truncated_distance = data["cleaved_distance"][row]

def runCmd(cmd):
    return run(session, cmd)

runCmd(f"open {pdb}")
runCmd("sel /" + chain)
runCmd("hide sel atoms")
runCmd("select ~sel")
runCmd("delete sel")
runCmd("delete solvent")
runCmd("delete ligand")
runCmd("delete ions")
runCmd("setattr model res_numbering uniprot")
runCmd(f"sel /{chain}:{original_position}")
runCmd("show sel atoms")

runCmd(f'open "/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/cleaved_N_Alphafold/Alphafold_cleaved_N_predictions/{gene_name_and_position}/{gene_name_and_position}_selected_prediction.pdb"')
runCmd(f"sel #2/A:{cleaved_position}")
runCmd("show sel atoms")
runCmd("sel @H*")
runCmd("hide sel")
runCmd("mm #2 to #1 showAlignment true")
runCmd("~select")

runCmd(f'label #1/{chain}:{original_position} text "{label}" color magenta')
# phi_symbol = "\u03D5"
# psi_symbol = "\u03C8"
runCmd(f"2dlabels text 'Original: relSESA={round(pdb_relSESA, 3)}, distance={round(pdb_distance, 3)}\u212B' xpos 0.02 ypos 0.94 bgColor white outline 1 margin 3")
runCmd(f"2dlabels text 'Cleaved: relSESA={round(truncated_relSESA, 3)}, distance={round(truncated_distance, 3)}\u212B' xpos 0.02 ypos 0.89 bgColor white outline 1 margin 3")

# runCmd("windowsize 1200 800")

runCmd("color #1 gold")
runCmd("color #1 byhetero")
runCmd("color #2 teal")
runCmd("color #2 byhetero")
runCmd('preset "overall look" "publication 2 (depth-cued)"')

print(gene_name_and_position)