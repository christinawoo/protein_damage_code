import pandas as pd
import json

df = pd.read_excel('/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/cleaved_n_docking/cleaved_n_docking_initial.xlsx')

df_dict = df.to_dict()

with open('/Users/matthewsu/Desktop/Chimera Computations/chimera/s24-extensions/cleaved_n_docking/preprocess_cleaved_n_docking.json', 'w', encoding='utf-8') as f:
    json.dump(df_dict, f, ensure_ascii=False, indent=4)
