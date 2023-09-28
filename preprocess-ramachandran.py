import pandas as pd
import json

df = pd.read_excel('pre-ramachandran.xlsx')

df_dict = df.to_dict()

# df_dict = df_dict[df_dict.columns.drop(list(df_dict.filter(regex='Unnamed:')))]

with open('new_ramachandran_preprocess.json', 'w', encoding='utf-8') as f:
    json.dump(df_dict, f, ensure_ascii=False, indent=4)

# f = open("preprocess.txt", "w")
# f.write(str(df_dict))
# f.close()
