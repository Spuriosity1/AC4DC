import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
from tbr import plotit, plt

plotit()
plt.savefig('/Users/alaric-mba/Desktop/tbr.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/qtbr_conservation.pgf')

