from scatter import res_to_q
import copy

print("====================")
print("Refreshed imaging params")
print("====================")

best_resolution = 2
worst_resolution = None
default_dict = dict( 

    #### Individual experiment arguments
    laser = dict(
        SPI = True,
        SPI_resolution = best_resolution,
        pixels_across = 100,  # for SPI, shld go on xfel params.
        random_orientation = True,  # infinite cryst sim only, NB: orientation is synced with undamaged crystal imaging  (TODO random orientation should not be separate from the other orient params...)
    ),
    ##### Crystal params
    crystal = dict(
        supercell_scale = 1,  ##3 # for SC: cell_scale^3 unit cells 
        num_supercells = 1,
        supercell_simulations = 1,        
        include_symmetries = True, ##True  # should unit cell contain symmetries?
        positional_stdv = 0, # Introduces disorder to positions. Note this is a deviation from the IDEAL structure, so is not a measure of similarity with undamaged and damaged structure but the ideal structure to recover and the dmaaged structure. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if quickly gauging serial crystallography R factor, as should somewhat average out.
        cell_packing = "SC",
        rocking_angle = 1.2,  # (approximating mosaicity, infinite crystal sim only)
        #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    ),
    #### XFEL params    
    experiment = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(worst_resolution),#None #angstrom
        q_cutoff = res_to_q(best_resolution),#2*np.pi/d
        t_fineness=25,   ##25
        #####crystal stuff
        max_miller_idx = None, # = m, where max momentum given by q with miller indices (m,m,m)
        ####SPI stuff
        num_rings = 20,
        pixels_per_ring = 25,
        # first image orientation (cardan angles, degrees) 
        SPI_x_rotation = 0,
        SPI_y_rotation = 0,
        SPI_z_rotation = 0,
        #crystallographic orientations (not consistent with SPI yet)
        num_orients_crys=20,
        orientation_axis_crys = None,#[1,1,0]   # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Incompatible with random_orientations   
    ),
)
asymmetric_unit_dict = copy.deepcopy(default_dict)
asymmetric_unit_dict["crystal"]["include_symmetries"] = False 



goldilocks_dict_unit = dict( 

    #### Individual experiment arguments
    laser = dict(
        SPI = True,
        SPI_resolution = best_resolution,
        pixels_across = 100,  # for SPI, shld go on xfel params.
        random_orientation = True,  # infinite cryst sim only, NB: orientation is synced with undamaged crystal imaging  (TODO random orientation should not be separate from the other orient params...)
    ),
    ##### Crystal params
    crystal = dict(
        supercell_scale = 1,  ##3 # for SC: cell_scale^3 unit cells 
        num_supercells = 1,
        supercell_simulations = 1,        
        include_symmetries = True, ##True  # should unit cell contain symmetries?
        positional_stdv = 0, # Introduces disorder to positions. Note this is a deviation from the IDEAL structure, so is not a measure of similarity with undamaged and damaged structure but the ideal structure to recover and the dmaaged structure. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if quickly gauging serial crystallography R factor, as should somewhat average out.
        cell_packing = "SC",
        rocking_angle = 1.2,  # (approximating mosaicity, infinite crystal sim only)
        #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    ),
    #### XFEL params    
    experiment = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(worst_resolution),#None #angstrom
        q_cutoff = res_to_q(best_resolution),#2*np.pi/d
        t_fineness=25,   ##25
        #####crystal stuff
        max_miller_idx = None, # = m, where max momentum given by q with miller indices (m,m,m)
        ####SPI stuff
        num_rings = 20,
        pixels_per_ring = 25,
        # first image orientation (cardan angles, degrees) 
        SPI_x_rotation = 0,
        SPI_y_rotation = 0,
        SPI_z_rotation = 0,
        #crystallographic orientations (not consistent with SPI yet)
        num_orients_crys=20,
        orientation_axis_crys = None,#[1,1,0]   # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Incompatible with random_orientations   
    ),
)
goldilocks_dict_3x3x3 = copy.deepcopy(goldilocks_dict_unit)
goldilocks_dict_3x3x3["crystal"]["supercell_scale"] = 3




# Water background
background_dict = copy.deepcopy(default_dict)
#background_dict["crystal"]["positional_stdv"] = 10  ### Idea: When we average over the distributions, this will form a background. 
background_dict["laser"]["SPI"] = True 
background_dict["crystal"]["positional_stdv"] = 0  ### Idea:  generate multiple solvents
background_dict["crystal"]["cell_scale"] = 1
background_dict["crystal"]["include_symmetries"] = False 
background_dict["experiment"]["t_fineness"] = 10  # We don't need to be as accurate with this.
##