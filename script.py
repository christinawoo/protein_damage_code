import json
import numpy as np
import os.path

# Working directory is desktop by default: https://www.cgl.ucsf.edu/chimerax/docs/user/commands/pwd.html

data = {}

with open(os.path.abspath("chimera/preprocess.json")) as json_file:
  data = json.load(json_file)

nan = data['relSESA']['0']
# np.isnan can be used to check! https://towardsdatascience.com/5-methods-to-check-for-nan-values-in-in-python-3f21ddd17eed
print(np.isnan(nan))
