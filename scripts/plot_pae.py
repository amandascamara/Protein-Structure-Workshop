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
pklfile = []
pklchain1 = []
pklchain2 = []
for l in range(1,len(lines)):
    pklfile.append(lines[l].rstrip('\n').split(' ')[0])
    pklchain1.append(lines[l].rstrip('\n').split(' ')[1])
    pklchain2.append(lines[l].rstrip('\n').split(' ')[2])

# reading the pdbfile 
parser = PDBParser()
structure = parser.get_structure("P", af2pdbfile)
uchains = []
clen = []
for chain in structure[0]:
    uchains.append(chain.get_id())
    clen.append(len(chain))
rchainid = []
for r in range(1,1+sum(clen)):
    for c in range(1,1+len(uchains)):
        if r>sum(clen[:c-1]) and r<=sum(clen[:c]): 
            rchainid.append(uchains[c-1])
# building maps
for f in range(len(pklfile)):
    print(pklfile[f])
    pkl = pickle.load(open(pklfile[f],'rb'))
    pae = pkl['predicted_aligned_error']
    paemap = np.zeros((len(structure[0][pklchain1[f]]),len(structure[0][pklchain2[f]])))
    a = -1
    for r1 in range(len(pae)):
        if rchainid[r1]==pklchain1[f]:
            a = a + 1
            b = -1
            for r2 in range(len(pae)):
                if rchainid[r2]==pklchain2[f]:
                    b = b + 1
                    paemap[a,b] = pae[a][b]
    
    plt.imshow(paemap, cmap="Greens_r", vmin=0, vmax=30)
    plt.colorbar(label="Expected position error (Ångströms)")
    folder = pklfile[f].rfind('/')+1
    folder2 = pklfile[f][:folder-1].rfind('/')+1
    plt.title(pklfile[f][folder2:])
    plt.xlabel("Scored residue of chain " + pklchain1[f],fontsize=16)
    plt.ylabel("Aligned residue of chain " + pklchain2[f],fontsize=16)
    figname = pklfile[f][folder2:folder-1] + '_' +  pklfile[f][folder:-4] + pklchain1[f] + pklchain2[f] + '_pae.png'
    plt.savefig(figname,dpi=300)
    plt.close()
