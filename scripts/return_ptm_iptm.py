from __future__ import print_function
import pickle
import sys

pkl = str(sys.argv[1])

p = pickle.load(open(pkl,'rb'))
print('The pTM is ', p['ptm'])
if len(p)==13: print('The interface pTM is ', p['iptm'])

