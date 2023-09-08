import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import sys

# color code
color_dict={"-":"0xffffff","H":"0x734F5A","B":"0x264653","E":"0x2A9D8F","G":"0xE9C46A","I":"0xE76F51","T":"0x941C2F","S":"0xC05761"}
# ss code
# - : irregular elements,
# H : alpha helix,
# B : isolated beta-bridge
# E : extended strand, participates in beta ladder
# G : 3/10 helix
# I : pi helix
# T : hydrogen bonded turn
# S : bend


# take chain info from each file
fileslist = str(sys.argv[1])

info = open(fileslist, "r")
print('Reading list of files and monomers.')
lines = info.readlines()
info.close()
filenames = []
for i in range(len(lines)):
    l = lines[i].rstrip('\n')
    if l[0]!='#':
        filenames.append(l.split(' ')[0])

print(filenames)

# read and write pymol command into pml file
for nf in range(len(filenames)):
    f = open(filenames[nf], "r")
    lines=f.readlines()
    f.close()
    nlines=len(lines)
    pmlfile = filenames[nf][:-4]+".pml"
    pml = open(pmlfile, 'w+')
    for l in lines:
        l=l.rstrip('\n')
        line = l.split(' ')[0]
        nres = l.split(' ')[1]
        chain = l.split(' ')[2]
        res = l.split(' ')[3]
        ss = l.split(' ')[4]
        pml.write("color %s, (resid %s and chain %s) \n" % (color_dict[ss],nres,chain))
    pml.close()
