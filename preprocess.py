import pandas as pd
import json

df = pd.read_excel('pre-ramachandran.xlsx')

df_dict = df.to_dict()

length = len(df_dict['gene_name_and_position'])

for i in range(length):
    index = str(i)
    df_dict["pdb_data_obj"][i] = json.loads(df_dict["pdb_data_obj"][i])

with open('preprocess.json', 'w', encoding='utf-8') as f:
    json.dump(df_dict, f, ensure_ascii=False, indent=4)

# f = open("preprocess.txt", "w")
# f.write(str(df_dict))
# f.close()
