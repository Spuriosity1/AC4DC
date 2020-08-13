import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BSpline
from math import log
import os.path as path
import sys
# import glob
import csv
import time
import subprocess

plt.style.use('seaborn')

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(9, 6))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Plotter:
    # Initialistation: Plotter(water)
    # expects there to be a control file named input/water.mol
    def __init__(self, mol):
        self.p = path.abspath(path.join(__file__ ,"../../"))
        # Inputs
        molfile = self.p+"/input/"+mol+".mol"
        self.mol = {'name': mol, 'file': molfile, 'mtime': path.getmtime(molfile)}
        # Stores the atomic input files read by ac4dc
        self.atomdict = {}

        # Outputs
        self.outDir = self.p + "/output/__Molecular/"+mol
        self.freeFile = self.outDir+"/freeDist.csv"
        self.intFile = self.outDir + "/intensity.csv"

        self.boundData={}
        self.freeData=None
        self.intensityData=None
        self.energyKnot=None
        self.timeData=None

        self.get_atoms()
        self.update()
        self.autorun=False

    # Reads the control file specified by self.mol['file']
    # and populates the atomdict data structure accordingly
    def get_atoms(self):
        self.atomlist = {}
        with open(self.mol['file'], 'r') as f:
            reading = False
            for line in f:
                if line.startswith("#ATOMS"):
                    reading=True
                    continue
                elif line.startswith("#"):
                    reading=False
                    continue
                if reading:
                    a = line.split(' ')[0].strip()
                    if len(a) != 0:
                        file = self.p + '/input/' + a + '.inp'
                        self.atomdict[a]={
                            'file': file,
                            'mtime': path.getmtime(file),
                            'out': self.outDir+"/dist_%s.csv"%a}

    def update_inputs(self):
        self.get_atoms()
        self.mol['mtime'] = path.getmtime(self.mol['file'])

    def rerun_ac4dc(self):
        cmd = self.p+'/bin/rates '+self.mol['file']
        print("Running: ", cmd)
        subprocess.run(cmd, shell=True, check=True)

    def check_current(self):
        # Pull most recent atom mod time
        self.update_inputs()
        # Outputs are made at roughly the same time: Take the oldest.
        out_time = path.getmtime(self.freeFile)
        for atomic in self.atomdict.values():
            out_time = min(out_time, path.getmtime(atomic['out']))

        if self.mol['mtime'] > out_time:
            print("Master file %s is newer than most recent run" % self.mol['file'])
            return False

        for atomic in self.atomdict.values():
            if atomic['mtime'] > out_time:
                print("Dependent input file %s is newer than most recent run" % atomic['file'])
                return False
        return True

    def parse_free_energy_spec(self):
        erow = []
        with open(self.freeFile) as f:
            r = csv.reader(f, delimiter=' ')
            for row in r:
                if row[0] != '#':
                    raise Exception('parser could not find energy grid specification, expected #  |')
                reading=False
                for entry in row:
                    # skip whitespace
                    if entry == '' or entry == '#':
                        continue
                    # it's gotta be a pipe or nothin'
                    if reading:
                        erow.append(entry)
                    elif entry == '|':
                        reading=True
                    else:
                        break
                if reading:
                    break
        return erow

    def update(self):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        self.intensityData = raw[:,1]
        self.timeData = raw[:, 0]
        self.energyKnot = np.array(self.parse_free_energy_spec(), dtype=np.float64)
        raw = np.genfromtxt(self.freeFile, comments='#', dtype=np.float64)
        self.freeData = raw[:,1:]
        for a in self.atomdict:
            raw = np.genfromtxt(self.atomdict[a]['out'], comments='#', dtype=np.float64)
            self.boundData[a] = raw[:, 1:]

    # makes a blank plot showing the intensity curve
    def setup_axes(self):
        fig.clf()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()
        ax2.plot(self.timeData, self.intensityData, lw = 2, c = 'black', ls = '--', alpha = 0.7)
        ax2.set_ylabel('Pulse fluence')
        ax1.set_xlabel("Time, fs")
        return (ax1, ax2)


    def plot_atom_raw(self, a):
        ax, pulseax = self.setup_axes()
        for i in range(self.boundData[a].shape[1]):
            ax.plot(self.timeData, self.boundData[a][:,i], label = 'Config '+str(i))
        ax.set_title("Charge state dynamics")
        ax.set_ylabel("Average charge")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

        plt.show()

    def plot_free(self):
        fig.clf()
        # Need to turn freeData (matrix of BSpline coeffs)
        # freeeData [t, c]
        # into matrix of values.
        num_E = 1000
        min = self.energyKnot[0]
        max = self.energyKnot[-1]
        Y = np.linspace(min, max, num_E)
        Z = np.zeros((num_E, self.timeData.shape[0]), dtype=np.float64)
        for t in range(self.freeData.shape[0]):
            inter = BSpline(self.energyKnot, self.freeData[t,:], 4)
            Z[:,t] = inter(Y)
        plt.contourf(self.timeData, Y, Z, cmap='RdGy')
        plt.title("Free electron energy distribution")
        plt.ylabel("Energy (eV)")
        plt.xlabel("Time, fs")
        plt.show()
        plt.colorbar()




pl = Plotter(sys.argv[1])
#pl.plot_free()
