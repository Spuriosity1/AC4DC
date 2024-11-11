#%%
'''
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================
(C) Spencer Passmore 2023
'''

import os
os.getcwd()
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
######

import os.path as path
from Bio.PDB.vectors import Vector as Bio_Vect
from Bio.PDB.vectors import homog_trans_mtx, set_homog_trans_mtx
from Bio.PDB.vectors import rotaxis2m
#from Bio.PDB.PDBParser import PDBParser
#from Bio.PDB.PDBIO import PDBIO
#from Bio.PDB.StructureBuilder import StructureBuilder
from xpdb import sloppyparser as xPDBParser
from xpdb import SloppyPDBIO as xPDBIO
from xpdb import SloppyStructureBuilder as xStructureBuilder  # Hack, enables atom counts over 10,000
from Bio.PDB.Atom import Atom as PDB_Atom
#from sympy.utilities.iterables import multiset_permutations
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import cos
from numpy import sin
from plotter_core import Plotter
from scipy.spatial.transform import Rotation as Rotation
from matplotlib.colors import to_rgb
import matplotlib as mpl
from matplotlib import cm
import copy
import pickle
import colorcet as cc; import cmasher as cmr
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.offline as pltly_offline
from IPython.display import display, HTML
from IPython import get_ipython


import pandas as pd
import csv

# Not implemented properly yet.

def read_reflection_file(fname):
    if fname[-4:] == ".rfl":
        fname = fname[:-4]
    directory = "reflections/"
    in_path = directory + fname + ".rfl"
    print("Reading reflection file",in_path)
    df = pd.read_csv(in_path,header=None,delim_whitespace=True,names=["h","k","l","I"])
    return df

fnames = ["hen_light_real_v2","hen_light_ideal_v2","hen_Gd_real_v2","hen_Gd_ideal_v2"]
num_gaps = 1000
binlims = np.linspace(0,1,num_gaps+1,endpoint=True)

cum_y = np.zeros(shape=(len(fnames),len(binlims))) 
for j,fname in enumerate(fnames):
    df_rfl = read_reflection_file(fname)
    df_rfl = df_rfl.sort_values(by='I',ascending=False)
    df_rfl = df_rfl.reset_index()
    df_rfl = df_rfl.drop(0) # drop (0,0,0)

    df_rfl["I"]/=df_rfl["I"].max()
    print(df_rfl)

    df_I = df_rfl["I"]

    for i, binlim in enumerate(binlims):
        cum_y[j][i] = df_I[df_I<binlim].count()
        #cum_y[j][i] = np.sum(df_I[df_I<binlim])

for i, dataset in enumerate(cum_y):
    plt.plot(binlims,dataset,label=fnames[i])
plt.show()

# %%
