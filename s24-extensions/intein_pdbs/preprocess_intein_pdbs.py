import pandas as pd
import json

df = pd.read_excel('/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/intein_pdbs/intein_pdbs_initial.xlsx')

df_dict = df.to_dict()

with open('/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/intein_pdbs/preprocess_intein_pdbs.json', 'w', encoding='utf-8') as f:
    json.dump(df_dict, f, ensure_ascii=False, indent=4)
