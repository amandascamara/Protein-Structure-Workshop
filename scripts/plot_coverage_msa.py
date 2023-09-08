### Script from https://elearning.bits.vib.be/courses/alphafold/ by Jasper Zuallaert (VIB-UGent), with the help of Alexander Botzki (VIB) and Kenneth Hoste (UGent).
### https://raw.githubusercontent.com/jasperzuallaert/VIBFold/main/visualize_alphafold_results.py

import glob
import math
import os
import numpy as np
from matplotlib import pyplot as plt
import sys
import pickle

def generate_output_images(feature_dict, out_dir):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    ##################################################################
    plt.figure(dpi=300)
    ##################################################################
    plt.title("Sequence coverage")
    plt.imshow(final,
	       interpolation='nearest', aspect='auto',
	       cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(f"{out_dir}/coverage.png")

    ##################################################################

input_dir = str(sys.argv[1])

feature_dict = pickle.load(open(f'{input_dir}/features.pkl','rb'))

generate_output_images(feature_dict, input_dir)

