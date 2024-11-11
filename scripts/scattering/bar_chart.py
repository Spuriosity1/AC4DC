#%%
import matplotlib.pyplot as plt
import numpy as np

target = ("Continuous protein", "Solvated crystal")
compositions = {
    'Light atoms': ([0.0177,0.0113],  [0.0250,0.0163]),  #[7.1 kev, 9 kev]
    'Lysozyme': ([0.0226,0.0150], [0.0266,0.0180]),
    'Lysozyme.Gd (10 mM Gd solvent)': ([0.0274,0.0306], [0.0312,0.0340]),
    'Lysozyme.Gd (water solvent)': ([0.0274,0.0306],[0.0294,0.0275])
}


x = np.arange(len(target))  # the label locations
width = 0.25  # the total width of the bars, before split by energy.
multiplier = 0

fig, ax = plt.subplots(layout='constrained')
cmap = plt.get_cmap("tab10")

handles = []
for attribute, energy_measurements in compositions.items():
    col = cmap(multiplier)
    for i, measurement in enumerate(energy_measurements):
        offset = width * multiplier + width/len(energy_measurements)*i
        rects = ax.bar(x + offset, measurement, width/len(energy_measurements), label=attribute,color=col)
        ax.bar_label(rects, padding=3)
    handles.append(rects)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Length (mm)')
ax.set_title('Penguin attributes by species')
ax.set_xticks(x + width, target)
ax.legend(handles=handles,loc='best', ncols=1)
ax.set_ylim(0, None)

plt.show()
# %%
import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    #thaumatin thing
    'axes.titlesize':10,     # fontsize of the axes title
    'axes.labelsize':10,    # fontsize of the x and y labels   
    'ytick.labelsize':10,
    'xtick.labelsize':10,
    'legend.fontsize':10,    
    'lines.linewidth':2,
})
import matplotlib.pyplot as plt
import numpy as np

COLUMNWIDTH = 3.4975
FIGWIDTH = COLUMNWIDTH*1.65#/2
FIGHEIGHT = FIGWIDTH*9/16#*9/16
DPI = 800        

# target = ("Continuous protein", "Solvated crystal")
# compositions = {
#     'Light atoms': ([0.0177,0.0113],  [0.0250,0.0163]),  #[7.1 kev, 9 kev]
#     'Lysozyme': ([0.0226,0.0150], [0.0266,0.0180]),
#     'Lysozyme.Gd (10 mM Gd solvent)': ([0.0274,0.0306], [0.0312,0.0340]),
#     'Lysozyme.Gd (water solvent)': ([0.0274,0.0306],[0.0294,0.0275])
# }


compositions = ("Light atoms", "Lsyozyme","Lys.Gd (10 mM Gd)","Lys.Gd (water)")
target = {
    "Continuous protein": ([0.0177,0.0113],[0.0226,0.0150],[0.0274,0.0306],[0.0274,0.0306],), #[7.1 kev, 9 kev]
    "Solvated crystal": ([0.0250,0.0163],[0.0266,0.0180],[0.0312,0.0340],[0.0294,0.0275]),
}
energy_dict = ["7.1","9.0"]


fig, ax = plt.subplots(layout='constrained')
cmap = plt.get_cmap("tab20")
colours = [cmap(4),cmap(5),cmap(0),cmap(1)]

num_splits = len(target)*len(energy_dict)
grouped_measurements =[None]*num_splits
labels = []
a = 0
# This is horrible and I apologise.
for attribute, structure in target.items():
    for c, composition in enumerate(structure):
        for e, energy_meas in enumerate(composition):
            if grouped_measurements[a*len(energy_dict) + e] is None:
                grouped_measurements[a*len(energy_dict) + e] = [] 
            grouped_measurements[a*len(energy_dict) + e].append(energy_meas)
            if c == 0:
                labels.append(attribute+", "+energy_dict[e]+" keV")
    a+=1
print(grouped_measurements)

x = np.arange(len(compositions))  # the label locations
width = 0.2  # the total width of the bars, before split by energy.
multiplier = 0
for attribute, measurement, i in zip(labels,grouped_measurements,range(len(compositions))):
    print(attribute)
    print(measurement)
    offset = width * multiplier
    col = colours[i]
    #col = [cmap(i) for i in range(len(compositions))]
    rects = ax.bar(x + offset, measurement,  width, label=attribute, color=col)
    #ax.bar_label(rects,padding=3)
            #ax.bar_label(rects, padding=3)
    
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('$R_{dmg}$')
ax.set_title('')
ax.set_xticks(x + width, compositions)
ax.legend(fancybox=True,ncol=2,loc='upper center',bbox_to_anchor=(0.5, 1.05),framealpha=1) #handletextpad=0.01,columnspacing=0.2,borderpad = 0.18

ax.set_ylim(0, 0.04)

plt.gcf().set_figwidth(FIGWIDTH)
plt.gcf().set_figheight(FIGHEIGHT)  
ax.yaxis.set_major_locator(plt.MaxNLocator(4))

#plt.show()
ax.ticklabel_format(style='sci', axis='y', scilimits=(-1,1),)
plt.savefig("barchart.png",dpi=DPI)
# %%
