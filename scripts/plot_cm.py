import numpy as np
import sys
from Bio.PDB import *
import matplotlib.pyplot as plt

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two) :
    """Returns a matrix of C-alpha distances between two chains"""
    answer = []
    r = -1
    for row, residue_one in enumerate(chain_one) :
        if residue_one.id[0]!="W" and residue_one.id[0][:2]!="H_":
            answer.append([])
            r = r + 1
            for col, residue_two in enumerate(chain_two) :
                if residue_two.id[0]!="W" and residue_two.id[0][:2]!="H_":
                    answer[r].append(calc_residue_dist(residue_one, residue_two))
    return answer

fileslist = str(sys.argv[1])

info = open(fileslist, "r")
print('Reading list of files and monomers.')
lines = info.readlines()
info.close()
filenames = []
chain1 = []
chain2 = []
for l in range(len(lines)):
    filenames.append(lines[l].rstrip('\n').split(' ')[0])
    chain1.append(lines[l].rstrip('\n').split(' ')[1])
    chain2.append(lines[l].rstrip('\n').split(' ')[2])

for nf in range(len(filenames)):
    if filenames[nf][-4:]=='.cif': parser = MMCIFParser()
    if filenames[nf][-4:]=='.pdb': parser = PDBParser()
    structure = parser.get_structure("P", filenames[nf])
    dist_matrix = calc_dist_matrix(structure[0][chain1[nf]], structure[0][chain2[nf]])
    contact_map = dist_matrix
    #contact_map = dist_matrix < 12.0

    plt.imshow(contact_map, cmap="Greens_r", vmin=0, vmax=12)
    plt.colorbar(label="Distance (Ångströms)")
    folder = filenames[nf].rfind('/')+1
    folder2 = filenames[nf][:folder-1].rfind('/')+1
    print(folder,folder2)
    plt.title(label = filenames[nf][:-4] + "_" + chain1[nf] + chain2[nf])
    plt.xlabel("Residue of chain " + chain1[nf],fontsize=16)
    plt.ylabel("Residue of chain " + chain2[nf],fontsize=16)
    if folder!=0: figname = filenames[nf][folder2:folder-1] + '_' +  filenames[nf][folder:-4] + '_' + chain1[nf] + chain2[nf] + '_cm.png'
    if folder==0: figname = filenames[nf][folder:-4] + '_' + chain1[nf] + chain2[nf] + '_cm.png'
    plt.savefig(figname,dpi=300)
    plt.close()
