# To extract 7-mer sequences around all N or Q in the human proteome

import numpy as np
import pandas as pd

swissprot = pd.read_csv(r'')
sequence = swissprot['Sequence'].tolist()
print(len(sequence)) 

length = len(sequence)
character = 'N' # to extract Q, simply replace 'N' to 'Q'
consensus = []
for i in range(length):
     for index in range(len(sequence[i])):
             if sequence[i][index]==character:
                    consensus.append(sequence[i][index-3:index+4])
print(len(consensus)) 


consensusallN = pd.DataFrame(consensus)
consensusallN.to_csv('', index=False)

Ndata = pd.read_csv(r'')
Nseq = Ndata['0'].tolist()

import random
randomN = random.sample(Nseq, 50000)
dfN = pd.DataFrame(randomN)
dfN.to_csv('', index=False)
