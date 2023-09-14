import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import csv

raw = np.genfromtxt(sys.argv[1],comments='#',dtype=np.float64)

erow=[]
with open(sys.argv[1]) as f:
    r = csv.reader(f, delimiter=' ')
    for row in r:
        if row[0]=='#' and row[1]=='|':
            erow = row
erow = np.array(erow[3:-1], dtype=np.float64)

t = raw[:,0]
f = raw[:,1:]
# f[t, e]

plt.clf()
plt.contourf(t, erow, f.T)
plt.ylabel('Energy (eV)')
plt.xlabel('time (fs)')
plt.colorbar()
plt.show()
