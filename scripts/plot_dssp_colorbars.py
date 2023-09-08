import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import sys

# color code
col_dict={0:"#ffffff",1:"#734F5A",2:"#264653",3:"#2A9D8F",4:"#E9C46A",5:"#E76F51",6:"#941C2F",7:"#C05761",8:"#c2c2c2"}
# ss code
# 0 : irregular elements,
# 1 : alpha helix,
# 2 : isolated beta-bridge
# 3 : extended strand, participates in beta ladder
# 4 : 3/10 helix
# 5 : pi helix
# 6 : hydrogen bonded turn
# 7 : bend
# 8 : end of chain


# take chain info from each file
fileslist = str(sys.argv[1])

info = open(fileslist, "r")
print('Reading list of files and monomers.')
lines = info.readlines()
info.close()
filenames = []
monomers = []
allmonomers = []
j=0
for i in range(len(lines)):
    l = lines[i].rstrip('\n')
    if l[0]!='#':
        dsspFileName = l.split(' ')[0]
        filenames.append(dsspFileName)
        nc = len(l.split(' '))
        print(i,nc)
        monomers.append([])
        for c in range(1,nc):
            monomers[j].append(l.split(' ')[c])
            allmonomers.append(l.split(' ')[c])
        j=j+1

ssdict=['-','H','B','E','G','I','T','S']
ssnumber=range(9)
umon=np.unique(allmonomers)
ssData = []
for i in range (len(umon)): ssData.append([])
foundMonomers=[-1]*len(umon)
print(foundMonomers)

print(filenames)
print(monomers)
print(umon)

# save filename for each monomer
fileLabels = []
for m in range(len(umon)):
    fileLabels.append([])
    for f in range(len(filenames)):
        for c in range(len(monomers[f])):
            if monomers[f][c]==umon[m]: fileLabels[m].append(filenames[f][:-9])
print(fileLabels)

# read and store ss data into lists
for nf in range(len(filenames)):
    f = open(filenames[nf], "r")
    lines=f.readlines()
    f.close()
    nlines=len(lines)
    line = []
    nres = []
    chain = []
    res = []
    ss = []
    for l in lines:
        l=l.rstrip('\n')
        line.append(l.split(' ')[0])
        nres.append(l.split(' ')[1])
        chain.append(l.split(' ')[2])
        res.append(l.split(' ')[3])
        ss.append(l.split(' ')[4])

#    print(np.unique(chain),np.unique(ss))

# update chain and ss values
    npchain=np.array(chain)
    _, idx = np.unique(npchain, return_index=True)
    unique_chains = npchain[np.sort(idx)]
    print(unique_chains)
    for nl in range(nlines):
        for s in range(8):
            if ss[nl]==ssdict[s]: ss[nl]=ssnumber[s]
    for c in range(len(unique_chains)):
        for m in range(len(umon)):
#            print(nf,c,m)
            if monomers[nf][c]==umon[m]:
                ssData[m].append([])
                foundMonomers[m] = foundMonomers[m] + 1
                for nl in range(nlines):
                    if chain[nl]==unique_chains[c]: ssData[m][foundMonomers[m]].append(ss[nl])
                if chain[nl]==unique_chains[c]: ssData[m][foundMonomers[m]].append(8)

print(ssData)

# check if all chains inside each monomers data have all the same length
for m in range(len(umon)):
    lengths=[]
    for f in range(foundMonomers[m]+1): lengths.append(len(ssData[m][f]))
    if not all(l==lengths[0] for l in lengths):
        print("Something strange with %s, some monomers have the same name but different sizes" % umon[m])
        print("Completing shorter monomers with null ss entry.")
        print(lengths)
    longer = np.amax(lengths)
    for f in range(foundMonomers[m]+1):
        while(len(ssData[m][f])<longer): ssData[m][f].append(8)

print(umon)
print(foundMonomers)
# plot
pallete = ListedColormap([col_dict[x] for x in col_dict.keys()])
print(pallete.N)
for m in range(len(umon)):
    fig, ax = plt.subplots()
    cax = ax.imshow(ssData[m], aspect=10, cmap=pallete)
    ax.set_yticks(range(foundMonomers[m]+1))
    ax.set_yticklabels(fileLabels[m])
    ax.tick_params(axis='y', which='major', labelsize=8)
    cbar = fig.colorbar(cax)
    cbar.ax.set_yticklabels(['Loop','Alpha helix','Beta bridge','Strand','Helix-3','Helix-5','Turn','Bend','null'])
    figname = umon[m] + '_dssp.png'
    plt.savefig(figname, dpi=300)

#    figure = plt.figure()
#    axes = figure.add_subplot(111)
#    axes.set_yticks(range(foundMonomers[m]+1))
#    axes.set_yticklabels(fileLabels[m])
#    axes.tick_params(axis='y', which='major', labelsize=8)
#    cm = axes.imshow(ssData[m], aspect=10, cmap=pallete)
#    cb = figure.colorbar(cm)
#    figname = umon[m] + '_dssp.png'
#    plt.savefig(figname, dpi=300)
    #plt.show()
