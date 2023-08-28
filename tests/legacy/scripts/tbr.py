import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 
import sys

plt.rc('axes', labelsize=12)    # fontsize of the x and y labels

rawQ = np.genfromtxt(sys.argv[1]+"_Q.csv")
rawGamma = np.genfromtxt(sys.argv[1]+"_gamma.csv")

energies = rawQ[:, 0]
scale = np.sqrt(np.outer(energies, energies))
Qdata = rawQ[:, 1:]/scale
Gdata = rawGamma[:, 1:]/scale



maxG = np.max(Gdata)
maxQ = np.max(Qdata)
bigg = max(maxQ,maxG)

#  colormap to use
cm = 'magma'

def plotit():    
    fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2, sharex=True, sharey=True, figsize = (2.5,2))
    bottom = 0
    
    im = ax1.pcolormesh(energies, energies, -Qdata, shading='nearest', vmin=bottom, vmax=bigg, cmap=cm,rasterized=True)
    ax1.set_title(r'Free')
    ax1.set_xlabel(r"$\epsilon_k$")
    ax1.yaxis.set_label_coords(-0.07, 0.07)
    ax1.set_ylabel(r"$\epsilon_l$",rotation=90)
    ax2.pcolormesh(energies, energies, Gdata, shading='nearest',vmin=bottom, vmax=bigg, cmap=cm,rasterized=True)
    # ax2.set_title(r'$\sum_{\eta} \Gamma^{TBR}_{kl}$',y=1.5)
    ax2.set_title(r'Bound')
    ax2.set_xlabel(r"$\epsilon_k$")
    ax2.xaxis.set_label_coords(0.1, -0.1)
    ax1.xaxis.set_label_coords(0.1, -0.1)
    
    # ax3.pcolormesh(energies, energies, np.abs(Gdata+ Qdata), vmin=bottom, vmax=bigg, cmap=cm)
    # ax3.set_title(r'$\left|\sum_{\eta} \Gamma^{TBR}_{kl} + \sum_j Q^{TBR}_{j, kl}\right|$',y=1.05)
    # ax4.pcolormesh(energies, energies, Gdata+0.5*Qdata, vmin=bottom, vmax=bigg, cmap=cm)
    # ax4.set_title(r'$\sum_{\eta} \Gamma^{TBR}_{kl} + \frac{1}{2}\sum_j Q^{TBR}_{j, kl}$',y=1.05)
    
    fig.subplots_adjust(left=0.13, right=0.95,top=0.88, bottom=0.4,wspace=0.1)
    cbar_ax = fig.add_axes([0.15, 0.25, 0.8, 0.02]) 
    
    cb = fig.colorbar(im, orientation="horizontal", cax=cbar_ax)
    cb.set_label(r'$R^{TBR}_{\xi, kl}$')
    


def plotsum():    
    fig, (ax) = plt.subplots(nrows=1,ncols=1, sharex=True, sharey=True, figsize = (3,3))
    im = ax.pcolormesh(energies, energies, (Gdata + Qdata), cmap = cm)
    ax.set_title('(G + Q)')
    fig.colorbar(im)

def plotdiff():    
    fig, (ax) = plt.subplots(nrows=1,ncols=1, sharex=True, sharey=True, figsize = (3,3))
    im = ax.pcolormesh(energies, energies, (Gdata + 0.5*Qdata), cmap = cm)
    ax.set_title('(G + 0.5Q)')
    fig.colorbar(im)

plotit()
# plotdiff()
# plotsum()
plt.show()
