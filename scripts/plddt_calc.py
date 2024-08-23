# For calculation of pLDDT score for the residue of interest in an AlphaFold predicted structure.

import re
import csv
import os

class HighPlddtResidue:
    def __init__(self, chain_id, residue_number, plddt_score):
        self.chain_id = chain_id
        self.residue_number = residue_number
        self.plddt_score = plddt_score

def extract_pLDDT_scores(filename):
    residues = []
    
    with open(filename, 'r') as file:
        for line in file:
            pattern = r'^[A-Z]\s[A-Z]{3}\s\d+\s+\d+\s+\d+\d+.*'
            if not re.match(pattern, line):
                continue
            
            parts = line.split()
            if len(parts) >= 5:
                chain_id = parts[0][0]
                residue_number = int(parts[2])
                plddt_score = float(parts[4])
                
                res = HighPlddtResidue(chain_id, residue_number, plddt_score)
                residues.append(res)
    return residues

def process_csv_and_cif(csv_file_path, cif_directory):
    results = []
    
    with open(csv_file_path, 'r') as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            uniprot_id = row['uniprot_id']
            aa_position = int(row['aa_position'])
            
            cif_file_name = f"AF-{uniprot_id}-F1-model_v4.cif"
            cif_file_path = os.path.join(cif_directory, cif_file_name)
            
            if not os.path.isfile(cif_file_path):
                print(f"CIF file not found for {uniprot_id}")
                continue
            
            residues = extract_pLDDT_scores(cif_file_path)
            
            for residue in residues:
                if residue.residue_number == aa_position:
                    is_high_confidence = "Yes" if residue.plddt_score >= 70 else "No"
                    results.append([uniprot_id, aa_position, residue.plddt_score, is_high_confidence])
                    break

    return results

def write_results_to_csv(results, output_csv_path):
    with open(output_csv_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['uniprot_id', 'aa_position', 'plddt_score', 'high_confidence'])
        writer.writerows(results)

if __name__ == "__main__":
    csv_file_path = 'your path.csv' 
    cif_directory = 'your directory' 
    output_csv_path = 'your output.csv'

    results = process_csv_and_cif(csv_file_path, cif_directory)
    write_results_to_csv(results, output_csv_path)
