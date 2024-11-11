import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plotter_core import Plotter
import sys
plt.rcParams["font.size"] = 16

label = sys.argv[1] +'_' + sys.argv[2]

pl = Plotter(sys.argv[1])

atomos = sys.argv[2]

minval = -10
maxval = 0

fig, ax = pl.plot_ffactor_get_R_AC4DC(atomos, 10, timespan=(minval,maxval), show_avg=False)
sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('plasma'), norm=plt.Normalize(vmin=minval, vmax=maxval))
 
cbaxes = fig.add_axes([0.5, 0.9, 0.35, 0.02])
a = plt.colorbar(sm, orientation='horizontal', cax = cbaxes)
cbaxes.set_xlabel('Time (fs)')

fig.set_size_inches(6,5)
fig.subplots_adjust(left=0.2, bottom=0.2 ,right=0.95, top=0.95)
fig.savefig('/Users/alaric-mba/Desktop/ffactor.png')
# fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_ffactor.pgf')
