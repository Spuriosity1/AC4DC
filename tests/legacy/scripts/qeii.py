import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

plt.style.use('seaborn-muted')

import sys
import csv


raw = np.genfromtxt(sys.argv[1],comments='#',dtype=np.float64)
fig, ax = plt.subplots(figsize = (2.5,2))
plt.subplots_adjust(left=0.25,bottom=0.24,right=0.95,top=0.95)
y1 = raw[:-2,2]
y2 = raw[:-2,3]
e = raw[:,1][:-2]
ax.plot(e[::3],(y2/np.sqrt(e))[::3],'.', label=r'Free')
ax.plot(e,y1/np.sqrt(e),'-', label=r'Bound')
ax.set_xlabel(r'Mean energy at index $k$ (Ha)')
ax.set_ylabel(r'$R_{\xi, k}$  ($a_0^{3}$Ha$\hbar^{-1}$)')
Y = np.abs(raw[:-2,2]-raw[:-2,3])
# Y = np.ma.masked_where(Y<=0,Y)

ax.tick_params(axis='y',which='minor')
# ax.set_ylim([1e-6,1e2])
# locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
# ax.yaxis.set_minor_locator(locmin)
# ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.plot(raw[:-2,1],Y,'k-',label='Difference')


ax.legend(loc='lower right')
# ax.set_yscale('log')

fig.savefig('/Users/alaric-mba/Desktop/eii.png')
fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/Qeii_conservation.pgf')
