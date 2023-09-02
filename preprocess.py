import pandas as pd
import json

df = pd.read_excel('short.xlsx')

df_dict = df.to_dict()

with open('preprocess.json', 'w', encoding='utf-8') as f:
    json.dump(df_dict, f, ensure_ascii=False, indent=4)

# f = open("preprocess.txt", "w")
# f.write(str(df_dict))
# f.close()

