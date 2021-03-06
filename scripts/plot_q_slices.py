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
from labellines import *

cmap = plt.get_cmap("tab10")

plt.rcParams["font.size"] = 16 # 8

label = sys.argv[1] +'_' + sys.argv[2]

# slices = [-2.5, 0, 2.5, 5]
slices = [-7.5, -5, -2.5, 0]
# slices = [0, 15, 23]
energies = [200, 500, 500, 1000]
colrs = [cmap(i) for i in range(4)]


pl = Plotter(sys.argv[1])

for (t, e, col ) in zip(slices, energies, colrs):
    lines = pl.plot_step(t, normed=True, color = col, lw=0.5)
    T = pl.plot_fit(t, e, normed=True, color=col, lw=0.5)

# pl.fig_steps.set_size_inches(3,2.5)
pl.fig_steps.set_size_inches(6,5)

pl.ax_steps.set_ylim([1e-4, 1])
pl.ax_steps.set_xlim([1,10000]) 

pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)
pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False

handles, labels = pl.ax_steps.get_legend_handles_labels()
order = [0,2,4,6,1,3,5,7]
# order = [0,2,4,1,3,5]
pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper left',ncol=2)
# labelLines(fitlines,align=True,xvals = [80,170,280,360])
plt.savefig('/Users/alaric-mba/Desktop/free_distribution_evolution.png')
# plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_slices.pgf')
