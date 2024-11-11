##### (C) Spencer Passmore 2023 
'''
File: interactive_handler.py
Purpose: Interface between user's desires and interactive_core.py.


'''
####
from types import SimpleNamespace 


#####

import os
import os.path as path
from interactive_core import InteractivePlotter#, fit_maxwell, maxwell
from _generate_plots import make_some_plots
import sys
import numpy as np
from QoL import set_highlighted_excepthook
sys.path.append(os.path.join(os.path.dirname(__file__), 'scattering'))
from core_functions import get_sim_params,get_sim_elements # Functions used for naming batches with the parameters we care about. 
import PIL
import io
from math import ceil

###### P_defaults
_P = SimpleNamespace()
_P.NORMALISE = False
_P.ELECTRON_DENSITY = False # if False, use energy density  
_P.ANIMATION = False # if true, generate animation rather than interactive figure (i.e. automatic slider movement) 
_P.ALSO_MAKE_PLOTS = False # Generate static plots for each simulation using _generate_plots.py
_P.SINGLE_FRAME = False # Save a png using plotly .
_P.NAMING_MODE = 0  # For legend. 0: full details of sim parameters + sim name | 1: elements in sim| 
_P.SCALE_DENSITY_BY_THOUSAND = False # Use cubic nm rather than cubic angstrom for measuring energy density  
_P.END_T = 9999  # Put at value to cutoff times early.
_P.POINTS = 70
if _P.ANIMATION:
    _P.POINTS = 200 # 200
_P.SINGLE_FRAME_DICT = dict(
    # In inches
    width = (7.24409436834*2/3),#3.49751*2*3/5,
    height = 7.24409436834/3, #3.49751*2/3
    line_width = 2,
    font_size = 10,
    superscript_font_size = 8,
    
    show_legend = False, 

    ylog = False, # < uses linear scale if False 
    xlog = False, # <

    #times = [-10],#[-10],#[0],#[-10,0], # fs.  Use multiple to plot traces from multiple times (colours will be same between frames however.)
    #y_range = [0,0.0014],#lin:[0,0.0014]/[0,0.009],#log: [-5,-1],      # < If axis is log scale, limits correspond to powers of 10.
    
    #times = [-10],
    #y_range = [0,0.0014], # -10 fs
    times = [0],
    y_range = [0,0.009], # 0 fs
    x_range = [None,7500], # None,         # < Use None for default. 
)
_P.INSET = False
_P.SUBPLOT_AX_FONT_COLOUR = "#495187"# "blue"
_P.SUBPLOT_AX_GRID_COLOUR = "#c4c4eb"
_P.INSET_DICT = dict(
    axes_kwargs = dict(
        xaxis2 =dict(
            domain=[0.1, 0.4],
            range = [0,120], #[0,50],
            anchor='y2',
            gridcolor=_P.SUBPLOT_AX_GRID_COLOUR, 
            zeroline = False,
            #zerolinecolor = _P.SUBPLOT_AX_GRID_COLOUR
        ),
        yaxis2=dict(
            domain=[0.15, 0.8],
            range = [0,None],
            anchor='x2',
            gridcolor=_P.SUBPLOT_AX_GRID_COLOUR, 
            zeroline = False,
            tickformat = ".1f",
            tick0=_P.SINGLE_FRAME_DICT['y_range'][0],
            dtick = (_P.SINGLE_FRAME_DICT['y_range'][1]- _P.SINGLE_FRAME_DICT['y_range'][0])/3*2*(1+(1000-1)*_P.SCALE_DENSITY_BY_THOUSAND),          
            #zerolinecolor = _P.SUBPLOT_AX_GRID_COLOUR
        ),
    ),
)
_P.TERMINAL_MODE = True
_P.MOLECULAR_PATH = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"

P_defaults = vars(_P)

######
def generate_graphs(param_dictionary,sys_argv=None):
    # Set up namespace by filling with default values
    if type(param_dictionary) is SimpleNamespace:
        param_dictionary = vars(param_dictionary)
    for key, default_value in P_defaults.items():
        if key not in param_dictionary.keys():
            param_dictionary[key] = default_value 
    assert P_defaults.keys() == param_dictionary.keys() 
    P = SimpleNamespace(**param_dictionary)
    assert P.ELECTRON_DENSITY is False, "Electron density not implemented yet, need to take in dynamic knots." #TODO
    if P.INSET:
        assert P.SINGLE_FRAME,"Inset only supported with single frame at present" 
    set_highlighted_excepthook()
    graph_folder = path.abspath(path.join(__file__ ,"../../output/_Graphs/")) + "/" 

    if sys_argv is None:
        sys_argv = sys.argv  
    if  len(sys_argv) < 2:
        print("Usage: python3 " +path.basename(__file__)+ " Carbon_1 <comparison_handle1...etc. (optional)>. Alternatively a folder name can be provided to generate interactives for all contained outputs.")
        exit()
   
    def make_outfolder(outdir,subfolders):
        for subdir in subfolders:
            try:
                os.makedirs(outdir+subdir,exist_ok=True)
            except OSError:
                print("Ignoring os.makedirs error.")
                pass
    
    def single_interactive(target_handles,data_parent_dir,outdir):
        print(target_handles,data_parent_dir,outdir)
        fname_out = target_handles[0]
        if len(sys_argv) > 2:
            fname_out += '_' + sys_argv[2][:15] + "..-"
        plot_title = fname_out.replace('_',' ')  
        int_subdir = "interactives/"; plot_subdir = "plots/"  # Separate out types of plots
        make_outfolder(outdir,[int_subdir,plot_subdir])
        if not P.SINGLE_FRAME:
            generate(P,target_handles,data_parent_dir,plot_title,fname_out + "_Int",outdir+int_subdir)  
        else:
            snapshot(P,target_handles,data_parent_dir,fname_out + "_frame",outdir+plot_subdir)  
        if P.ALSO_MAKE_PLOTS:
            for i, target in enumerate(target_handles):
                if len(sys_argv) > 2:
                    fname_out = target
                make_some_plots(target,data_parent_dir,fname_out+"_Plt",outdir+plot_subdir,True,True,False,False)   #TODO bound plots not working for multiple species, they are joined with multiple species.          

    given_path = ""
    is_parent_dir = False
    if len(sys_argv) == 2:
        is_parent_dir = True
        given_path = P.MOLECULAR_PATH + sys_argv[1] + "/"
        if path.isfile(given_path):
            raise Exception("Please pass simulation output folder name or folder containing simulation output folders")
        for f in os.listdir(given_path):
            if path.isfile(given_path + "/"+ f):
                is_parent_dir = False
                break
    if is_parent_dir:
        # Generate interactive for each output contained in the parent directory
        outdir = graph_folder + sys_argv[1] + "/"
        for target_handle in os.listdir(given_path):
            single_interactive([target_handle],given_path,outdir)
        # ipl = multi_interactive(given_path)
        #save_interactive(ipl,outdir,sys_argv[1] + "_interactive")
        print("Interactives successfully saved to",outdir)
    else: 
        outdir = graph_folder
        target_handles = []  
        # Store folder names from commandline args without .mol extension. 
        for i, target in enumerate(sys_argv):
            if i == 0: continue #####
            if target.endswith(".mol"):
                target = target[0:-4]         
            target_handles.append(target)        
        single_interactive(target_handles,P.MOLECULAR_PATH,outdir)
    print("Done! Remember to drink water!")


def set_up_interactive_axes(P,ipl,plot_title,font_size=35,superscript_font_size=30):
    '''
    superscript_font_size is for whatever plotly does when using superscript tick labels, font is made bigger and I assume it bases it off the superscript or something.
    '''

    # Y scaling
    lin_ymax = 0.03
    lin_ymin = -0.005
    log_ymax = 1
    log_ymin = 1e-4
    if P.ANIMATION:
        log_ymax = 0.5
        log_ymin = 5e-5
    xmin,xmax = ipl.get_E_lims()
    xmin = max(1,xmin)
    xmax/=2.4 # Cut off tail

    if P.SCALE_DENSITY_BY_THOUSAND:
        P.SINGLE_FRAME_DICT["y_range"] =[v*1000 for v in P.SINGLE_FRAME_DICT["y_range"]] 
    # TODO pop in units
    ylabel = "Electron energy density" 
    if P.ANIMATION or P.SINGLE_FRAME:
        ylabel = "Energy density"
    if P.ELECTRON_DENSITY:
    # Latex in plotly is broken. 
        ylabel = "Electron density"
        log_ymin = 2e-7; log_ymax = 1e-2
        if not P.NORMALISE and not P.ANIMATION:
            #ylabel = "$\\text{"+ylabel+" (} \\AA^{-3}"
            if P.SCALE_DENSITY_BY_THOUSAND:
                ylabel +=" (eV nm^{-3})" 
            else: 
                ylabel +=" (ev ang^{-3})"
               # ylabel += " \\texttimes 1000 "
           # ylabel += "$)"
    else:
        if not P.NORMALISE:
            #ylabel = "$\\text{"+ylabel+" (} "
            if P.SCALE_DENSITY_BY_THOUSAND:
                ylabel += " (eV nm^{-3})" 
                #ylabel+="keV\\cdot "
            #else:
                #ylabel +="eV\\cdot "
            #ylabel += "\\text{r{A}}^{-3}" 
            #ylabel+=")$"
        
    # Axis parameters, these define axis properties depending on what scale is used.
    xlabel = 'Energy (eV)'
    #ylabel = '$f(\\epsilon) \\Delta \\epsilon$' 
    if P.NORMALISE: 
        ylabel += " (arb. units)"
    

    if P.ELECTRON_DENSITY:
        ylabel = '$f(\\epsilon)'
    #x_lin_args = {'title': {"text": xlabel + " - lin scale", "font":{"size": 30,"family": "times new roman"}}, 'tickfont': {"size": 20}, 'type' : "linear", "range" : [xmin,xmax]}
    #y_lin_args = {'title': {"text": ylabel + " - lin scale", "font":{"size": 30,"family": "times new roman"}}, 'tickfont': {"size": 20}, 'type' : "linear", "range" : [lin_ymin,lin_ymax]}
    #x_log_args = {'title': {"text": xlabel + " - log scale", "font":{"size": 30,"family": "times new roman"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
    #y_log_args = {'exponentformat':'e','tick0': [100,10,1,0.1,0.01,0.001,0.0001,0.00001],'title': {"text": ylabel + " - log scale", "font":{"size": 30,"family": "times new roman"}}, 'tickfont': {"size": 20}, 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}
    
    x_lin_args = {'tick0': None,'dtick':None,'title': {"text": xlabel, "font":{"size": font_size,"family": "times new roman"}}, 'tickfont': {"size": font_size}, 'type' : "linear", "range" : [xmin,xmax]}
    y_lin_args = {'tick0': None,'dtick':None,'title': {"text": ylabel, "font":{"size": font_size,"family": "times new roman"}}, 'tickfont': {"size": font_size}, 'type' : "linear", "range" : [lin_ymin,lin_ymax]}
    x_log_args = {'tick0': [2,1,0,-1,-2,-3,-4,-5],'dtick':"1",'exponentformat':'power','showexponent':'all','title': {"text": xlabel, "font":{"size": font_size,"family": "times new roman"}}, 'tickfont': {"size": superscript_font_size,"family": "times new roman"}, 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
    y_log_args = {'tick0': [2,1,0,-1,-2,-3,-4,-5],'dtick':"1",'exponentformat':'power','showexponent':'all','title': {"text": ylabel, "font":{"size": font_size,"family": "times new roman"}}, 'tickfont': {"size": superscript_font_size,"family": "times new roman"}, 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}


    #for k,lab in zip(["yaxis2","xaxis2"],[ylabel,xlabel,]):
    #    P.INSET_DICT["axes_kwargs"][k]["text"] = lab
    P.INSET_DICT["axes_kwargs"]["xaxis2"]["title"] ={"text": xlabel, "standoff":0, "font":{"color":P.SUBPLOT_AX_FONT_COLOUR, "size": font_size,"family": "times new roman"}}

    if P.ANIMATION:
        x_log_args = {'tick0': [2,1,0,-1,-2,-3,-4,-5],'dtick':"1",'exponentformat':'power','showexponent':'all','title': {"text": xlabel, "font":{"size": font_size*2-5,"family": "times new roman"}}, 'tickfont': {"size": superscript_font_size*2-5,"family": "times new roman"}, 'type' : "log", "range" : [np.log10(xmin),np.log10(xmax)]}
        y_log_args = {'tick0': [2,1,0,-1,-2,-3,-4,-5],'dtick':"1",'exponentformat':'power','showexponent':'all','title': {"text": ylabel, "font":{"size": font_size*2-5,"family": "times new roman"}}, 'tickfont': {"size": superscript_font_size*2-5,"family": "times new roman"}, 'type' : "log", "range" : [np.log10(log_ymin),np.log10(log_ymax)]}
    if P.SINGLE_FRAME:
        y_lin_args['tick0']=P.SINGLE_FRAME_DICT['y_range'][0] 
        y_lin_args['dtick'] = (P.SINGLE_FRAME_DICT['y_range'][1]- P.SINGLE_FRAME_DICT['y_range'][0])/3
        y_lin_args['exponentformat']='power' # does nothing?
        y_lin_args['tickformat'] = ".1f"
    # Initialises graph object.
    ipl.initialise_figure(plot_title, x_log_args,y_log_args) 
    return (x_log_args,x_lin_args,y_log_args,y_lin_args)

def get_legend_labels(P,target_handles,sim_data_parent_dir,sim_input_dir=None):
    if sim_input_dir is None:
        sim_input_dir = path.abspath(path.join(__file__,"../../input/")) + "/"        
    custom_names = []
    legend_title=None
    if P.NAMING_MODE == 1:
        legend_title = "Elements present"        
    for handle in target_handles:
        if P.NAMING_MODE == 0:
            s = get_sim_params(handle,sim_input_dir,sim_data_parent_dir)
            sim_params = [s[0]["energy"],s[0]["width"],s[0]["fluence"]]
            unit_list = s[1]
            param_list = s[2]
            name = ["  " + param_list[i][0] + ": " +str(p) + unit_list[i]  for i, p in enumerate(sim_params)]
            name = ''.join(name) + " [" + handle + "]"
            name = name[2:]  
        elif P.NAMING_MODE == 1:
            elements = list(get_sim_elements(handle,sim_input_dir,sim_data_parent_dir))
            # Remove suffixes
            for e,elem in enumerate(elements):
                elements[e] = elem.split("_")[0]
            # Combine light atoms 
            light_atoms = ['C', 'N', 'O']
            if all(X in elements for X in light_atoms):
                for X in light_atoms:
                    elements.remove(X)
                elements.insert(0,"CNO")
            name = ','.join(elements)
        custom_names.append(name)
    return custom_names, legend_title

def generate(P,target_handles,sim_data_parent_dir,plot_title,fname_out,outdir):
    '''
    Arguments:
    target_handles: The name of the folder containing the simulation's data (the csv files). 
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''
    print("[ Interactive ] Creating electron density interactive figure")
    custom_names,legend_title = get_legend_labels(P,target_handles,sim_data_parent_dir)
    presentation_mode = None # Use default
    if P.ANIMATION:
        presentation_mode = True
    # Initialises plotter object with data from files.
    ipl = InteractivePlotter(target_handles,sim_data_parent_dir, max_final_t=P.END_T,max_points=P.POINTS,custom_names=custom_names,use_electron_density = P.ELECTRON_DENSITY,presentation_mode=presentation_mode,inset=P.INSET)  
    scale_button_args = set_up_interactive_axes(P,ipl,plot_title)
    # The meat of the plotting. Plot line for each point in time
    line_1 = {'width': 6,"dash": '10,1'}
    line_2 = {'width': 6,}
    line_kwargs = [line_1,line_2]
    if len(target_handles) != 2:
        line_kwargs = [line_2 for l in range(len(target_handles))]
    ipl.plot_traces(normed=P.NORMALISE, line_kwargs=line_kwargs)
    # Add widgets
    ipl.add_time_slider(one_slider=P.ANIMATION)  
    if not P.ANIMATION:
        ipl.add_scale_button(*scale_button_args)                 
    #-----Save-----#
    if P.ANIMATION:
        save_animation(ipl,outdir,fname_out,max_num_frames=999)
    else:
        save_interactive(ipl,outdir,fname_out)

def snapshot(P,target_handles,sim_data_parent_dir,fname_out,outdir):
    '''
    Arguments:
    target_handles: The name of the folder containing the simulation's data (the csv files). 
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''
    print("[ Interactive ] Creating snapshot of electron density interactive figure")
    custom_names,legend_title = get_legend_labels(P,target_handles,sim_data_parent_dir)
    font_size = P.SINGLE_FRAME_DICT["font_size"]
    supfont_size = P.SINGLE_FRAME_DICT["superscript_font_size"]

    
    ipl = InteractivePlotter(target_handles,sim_data_parent_dir, font_size=font_size, max_final_t=-99999,max_points=1,legend_title=legend_title,times_in_legend=False,custom_names=custom_names,use_electron_density = P.ELECTRON_DENSITY,inset=P.INSET)
    # Axes
    xlog_args,xlin_args,ylog_args,ylin_args = set_up_interactive_axes(P,ipl,None,font_size=font_size,superscript_font_size=supfont_size)
    x_args = xlin_args
    y_args = ylin_args 
    if P.SINGLE_FRAME_DICT["xlog"]:
        x_args = xlog_args          
    if P.SINGLE_FRAME_DICT["ylog"]:
        y_args = ylog_args 

    # Adjust layout to suit a static figure.

    for i, key in enumerate(["x_range","y_range"]):
        for j in range(2):
            val = P.SINGLE_FRAME_DICT[key][j]
            if val is not None:
                if i == 0:
                    x_args["range"][j] = val
                if i == 1:
                    y_args["range"][j] = val
            


    ipl.fig.update_layout(
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0,
            pad=0
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor = "left",
            x = 0.05,
            bgcolor = 'rgba(255,255,255,0.7)',
        ),          
        autosize=False,
        width=P.SINGLE_FRAME_DICT["width"]*72.27,
        height=P.SINGLE_FRAME_DICT["height"]*72.27,
    )   
    ipl.fig.update_xaxes(x_args)
    ipl.fig.update_yaxes(y_args)    

    for t in P.SINGLE_FRAME_DICT["times"]:
        # Plot a trace for each simulation for the snapshot in time. 
        ipl.input_data_args["max_final_t"] = t
        ipl.initialise_data()
        if P.SCALE_DENSITY_BY_THOUSAND:
            for i in range(len(ipl.target_data)):
                ipl.target_data[i].freeData*=1000
        line = {'width': P.SINGLE_FRAME_DICT["line_width"]}
        ipl.plot_traces(normed=P.NORMALISE, line_kwargs=[line]*len(target_handles))   
  
    ipl.fig.update_layout(
        showlegend=P.SINGLE_FRAME_DICT["show_legend"]
    )

    P.INSET_DICT["axes_kwargs"]["yaxis2"].update(
        #1:1 yaxis scale
        #range = [P.INSET_DICT["axes_kwargs"]["yaxis2"]["domain"][0]*P.SINGLE_FRAME_DICT["y_range"][1], P.INSET_DICT["axes_kwargs"]["yaxis2"]["domain"][1]*P.SINGLE_FRAME_DICT["y_range"][1]],
        #showticklabels = False,

        ####
        tickfont = y_args["tickfont"],
        side= 'right',        
    )
    P.INSET_DICT["axes_kwargs"]["xaxis2"].update(
        tickfont = x_args["tickfont"],
        side= 'top',
    )
    for k in ["yaxis2","xaxis2"]:
        P.INSET_DICT["axes_kwargs"][k]["tickfont"]["color"] = P.SUBPLOT_AX_FONT_COLOUR


    if P.INSET:
        ipl.fig.update_layout(**P.INSET_DICT["axes_kwargs"])
        ipl.fig.update_layout(**P.INSET_DICT["axes_kwargs"])
        if P.SINGLE_FRAME_DICT["show_legend"]:
            ipl.fig.update_layout(
                legend=dict(
                    xanchor = "left",
                    x = P.INSET_DICT["axes_kwargs"]["xaxis2"]["domain"][1]+0.05,
                ), 
            )
    
    
    if len(P.SINGLE_FRAME_DICT["times"]) == 1:
        fname_out +=  "_" + str(P.SINGLE_FRAME_DICT["times"][0])+"fs"
    save_snapshot(ipl,outdir,fname_out)

def save_interactive(ipl,outdir,fname_out):
    extension = ".html"
    file_path =  outdir + fname_out + extension
    ipl.fig.write_html(file_path)

    print("Done!") 

def save_snapshot(ipl,outdir,fname_out):
    import time
    import plotly.express as px
    file_format = "svg"
    file_path = outdir + fname_out + "."+file_format
    # Workaround for loading mathjax graphical bug
    ipl.fig.write_image(file_path,format=file_format)  
    time.sleep(2)

    ipl.fig.write_image(file_path,format=file_format)

    print("done!")

def save_animation(ipl,outdir,fname_out,max_num_frames=10):
    # https://stackoverflow.com/questions/55460434/how-to-export-save-an-animated-bubble-chart-made-with-plotly
    extension = ".gif"
    file_path =  outdir + fname_out + extension

    # generate images for each step in animation
    frames = []

    steps = ipl.steps_groups[-1] 
    #for i, step in enumerate(steps[::ceil(len(steps)/num_frames)]):
    num_frames = min(max_num_frames,len(steps))  
    #for i in range(len(steps)):  
    for i in range(len(steps)-35): #TODO remove, use above line instead.  
        if i%ceil(len(steps)/num_frames) != 0: continue # Iterate through every nth step to get our frame
        print("step=",i)
        # # set main traces to appropriate traces within plotly frame
        # ipl.fig.update(data=fr.data)
        # move the slider containing all simulation traces to correct place
        ipl.fig.layout.sliders[-1].update(active=i)  
        # For some reason this isn't making the trace update...
        #ipl.fig.update_layout(style = {"visible": ipl.steps_groups[-1][i]["args"][0]["visible"]}) 

        # fig.for_each_trace(
        #     lambda trace: trace.update(visible=True) if trace.name == "linear" else (),
        # )        
        visible = ipl.steps_groups[-1][i]["args"][0]["visible"]
        ipl.fig.for_each_trace(
            lambda trace: trace.update(visible=False),
        )
        for elem in np.where(visible)[0]:  # Get indices of traces that are to be visible at slider step.
            ipl.fig.data[elem].visible = True
        # generate image of current state
        frames.append(PIL.Image.open(io.BytesIO(ipl.fig.to_image(format="png", width=2560, height=1440, scale=1))))

    # append duplicated last image more times, to keep animation stop at last status
    for i in range(ceil(num_frames/5)):
        frames.append(frames[-1])

    # # Duplicate initial images for some Sanders-patented tension. 
    # for i in range(int(num_frames/10)):
    #     frames.insert(0,frames[0])

    # create animated GIF
    duration = 10 # duration in seconds
    frames[0].save(
            file_path,
            save_all=True,
            append_images=frames[1:],
            optimize=False,
            duration=duration*1000/num_frames,
            loop=0,
        )

if __name__ == "__main__":
    generate_graphs(P_defaults)


############### Test defaults
#python3.9 scripts/generate_interactive.py my_mol_1