import numpy as np
import sys
from Bio.PDB import *
import matplotlib.pyplot as plt

# take chain info from each file
fileslist = str(sys.argv[1])

info = open(fileslist, "r")
print('Reading list of files and monomers.')
lines = info.readlines()
info.close()
filenames = []
chains = []
for i in range(len(lines)):
    l = lines[i].rstrip('\n')
    if l[0]!='#':
        filenames.append(l.split(' ')[0])
        chains.append(l.split(' ')[1])
    print(filenames[i],chains[i])

for nf in range(len(filenames)):
    bfactors = []
    resid = []
    if filenames[nf][-4:]=='.cif': parser = MMCIFParser()
    if filenames[nf][-4:]=='.pdb': parser = PDBParser()
    structure = parser.get_structure("P", filenames[nf])
    for residue in structure[0][chains[nf]].get_residues():
        if residue.id[0]!="W" and residue.id[0][:2]!="H_":
            bfactors.append(residue["CA"].get_bfactor())
            resid.append(residue.id[1])
    label = filenames[nf][:-4] + "_" + chains[nf]
    plt.plot(resid,bfactors,label=label,linewidth=2)

plt.xlabel("Residue",fontsize=16)
plt.ylabel("B factor",fontsize=16)
plt.legend()
plt.savefig('Bfactors.png',dpi=300)
