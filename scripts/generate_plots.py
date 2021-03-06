import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plot_molecular_charge import Plotter
import sys

label = sys.argv[1] +'_' + sys.argv[2]

pl = Plotter(sys.argv[1])

pl.plot_tot_charge(every=10)
pl.fig.set_size_inches(3,2.5)
plt.savefig('/Users/alaric-mba/Desktop/charge_conservation.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_qcons.pgf')

pl.fig.set_size_inches(3,2.5)
pl.plot_free(log=True, min=1e-7, every=5)
pl.fig_free.subplots_adjust(left=0.15)
plt.savefig('/Users/alaric-mba/Desktop/free.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_free.pgf')

pl.plot_all_charges()
pl.fig.set_size_inches(3,2.5)
pl.fig.subplots_adjust(left=0.2,bottom=0.18,top=0.95)
plt.savefig('/Users/alaric-mba/Desktop/dynamics.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_bound.pgf')