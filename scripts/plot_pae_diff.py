import numpy as np
import sys
from Bio.PDB import *
import matplotlib.pyplot as plt
import pickle

def get_pae_plddt(model_names):
    out = {}
    for i,name in enumerate(model_names):
        d = pickle.load(open(name,'rb'))
        basename = os.path.basename(name)
        basename = basename[basename.index('model'):]
        out[f'{basename}'] = {'plddt': d['plddt'], 'pae':d['predicted_aligned_error']}
    return out

fileslist = str(sys.argv[1])

info = open(fileslist, "r")
print('Reading list of files and monomers.')
lines = info.readlines()
info.close()
af2pdbfile = lines[0].rstrip('\n')
pklfile1 = lines[1].rstrip('\n').split(' ')[0]
pkl1chain1 = lines[1].rstrip('\n').split(' ')[1]
pkl1chain2 = lines[1].rstrip('\n').split(' ')[2]
pklfile2 = lines[2].rstrip('\n').split(' ')[0]
pkl2chain1 = lines[2].rstrip('\n').split(' ')[1]
pkl2chain2 = lines[2].rstrip('\n').split(' ')[2]
pkl1 = pickle.load(open(pklfile1,'rb'))
pae1 = pkl1['predicted_aligned_error']
pkl2 = pickle.load(open(pklfile2,'rb'))
pae2 = pkl2['predicted_aligned_error']

if len(pae1)!=len(pae2): print("Something wrong! The two files don't contain the same number of residues.")

parser = PDBParser()
structure = parser.get_structure("P", af2pdbfile)
uchains = []
clen = []
for chain in structure[0]:
    uchains.append(chain.get_id())
    clen.append(len(chain))
 

rchainid = []
for r in range(1,1+len(pae1)):
    for c in range(1,1+len(uchains)):
        if r>sum(clen[:c-1]) and r<=sum(clen[:c]): 
            rchainid.append(uchains[c-1])

paemap = np.zeros((len(structure[0][pkl1chain1]),len(structure[0][pkl1chain2])))
a = -1
for r1 in range(len(pae1)):
    if rchainid[r1]==pkl1chain1:
        a = a + 1
        b = -1
        for r2 in range(len(pae1)):
            if rchainid[r2]==pkl1chain2:
                b = b + 1
                paemap[a,b] = pae1[a][b] - pae2[a][b]

plt.imshow(paemap, cmap="bwr", vmin=-5, vmax=5)
plt.colorbar(label="Expected position error (Ångströms) \n difference between model1 and model2")
plt.xlabel("Scored residue of chain " + pkl1chain1,fontsize=16)
plt.ylabel("Aligned residue of chain " + pkl1chain2,fontsize=16)
plt.savefig("pae_diff.png",dpi=300)

