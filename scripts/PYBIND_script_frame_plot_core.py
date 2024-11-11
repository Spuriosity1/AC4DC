'''===========================================================================
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
===========================================================================*/
// (C) Spencer Passmore 2023

[!!] WARNING: integrated with plasma simulation, edit at own risk! [!!]
'''

import numpy as np
# from scipy.interpolate import BSpline
from math import log
import os.path as path
import os
import sys
# import glob
import csv
import subprocess
import re
import random
from scipy.optimize import curve_fit
import plotly.graph_objects as go


class Plotter:
    # Example initialisation: Plotter(water)
    # --> expects there to be a control file named water.mol within AC4DC/input/ or a subdirectory.
    def initialise_figure(self, target, x_args={}, y_args={}):
        self.x_args = x_args
        self.y_args = y_args
        self.fig = go.Figure()
        
        self.fig.update_layout(
            title= target + " - Free-electron distribution",  # Attention: title overwritten by add_time_slider()
            showlegend=False,
            font=dict(
                family="Courier New, monospace",
                size=20,
                # color='rgba(0,0,0,0)' # hide tick numbers? Nope.
            ),
            paper_bgcolor= '#F7CAC9', #'#f7dac9' '#F7CAC9'  '#FFD580' "#92A8D1"  lgrey = '#bbbfbf',
            width=1920,
            height=1080            
        )
        self.fig.update_xaxes(x_args)
        self.fig.update_yaxes(y_args)  

        self.density = []
        self.energyKnot = []

        # initialise trace
        self.fig.add_trace(go.Scatter(x=[0,1],y=[0,1],name="density"))

    def update_data(self,energies,densities):
        self.density = np.array(densities)
        self.energyKnot = np.array(energies)         
    # Plots a single point in time.
    def plot_step(self, normed=True):        
        # norm = np.sum(self.freeData[n,:])
        dens = self.density
        E = self.energyKnot

        if normed:
            tot = self.get_density()
            dens /= tot
            dens /= 4*3.14

        self.fig.update_traces(x=E,y=dens*E)        
    def plot_the_knot(self,knot_to_plot):
        #self.fig.add_scattergl(x=knot_to_plot,y = [0.9]* len(knot_to_plot))    
        self.fig.add_scatter(x=knot_to_plot,y = [10**(self.y_args["range"][0])*1.1]* len(knot_to_plot),mode="markers",marker = dict(color='#e66000',size=8),name="knot")   
        self.fig.update_layout(
            showlegend=True,
            legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor = '#F5F5F5'
        ))         

    def plot_fit(self, fitE, normed=True, **kwargs):
        fit = self.energyKnot.searchsorted(fitE)
        dens = self.density
        if normed:
            tot = self.get_density()
            dens /= tot
            dens/=4*3.14

        Xdata = self.energyKnot[:fit]
        Ydata = dens[:fit]
        T, n = fit_maxwell(Xdata, Ydata)
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, T, n)*self.energyKnot,
            '--',label='%3.1f eV' % T, **kwargs)

    def plot_maxwell(self, kT, n, **kwargs):
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, kT, n)*self.energyKnot,
            '--',label='%3.1f eV' % kT, **kwargs)

    def get_density(self):
        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])
        de = de [1:] - de[:-1]
        return np.dot(self.density, de)

        

def fit_maxwell(X, Y):
    guess = [200, 12]
    # popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess, sigma=1/(X+10))
    popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess)
    return popt

def maxwell(e, kT, n):
    if kT < 0:
        return 0 # Dirty silencing of fitting error - note we get negative values from unphysical oscillations, so this increases the average value around this point. -S.P.
    return n * np.sqrt(e/(np.pi*kT**3)) * np.exp(-e/kT)

def lnmaxwell(e, kT, n):
    return np.log(n) + 0.5*np.log(e/np.pi*kT**3) - e /kT

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
