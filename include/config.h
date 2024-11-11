/**
 * @file config.h
 * @brief 
 */

#pragma once

const int GLOBAL_BSPLINE_ORDER = 3;  // 1 = rectangles, 2=linear, 3=quadratic  ATTENTION I have not yet gone through and ensured that there are no local redefinitions of the order.

#ifdef DEBUG
#define OUTPUT_TBR_TO_CERR
#define OUTPUT_DQDT_TO_CERR
#define OUTPUT_DFDT_TO_CERR
#endif


//-Disable plasma processes  
//#define NO_TBR    //Three body recombination
//#define NO_EE   // Electron-electron scattering. This seems to break the dynamic grid late in the simulation depending on pulse parameters.
//#define NO_EII    // Electron impact ionisation0
//#define NO_PLASMA // Disables all of the above. (i.e. primary - Photo,Auger,Fluoro, - ionisation only)

/// Disable features
//#define NO_PLOTTING // Turns off live saves of the free-electron energy distribution to _live_plot.png. Disables use of python 
//#define NO_BACKUP_SAVING // Disables the hourly saves of the data to  
//#define TRACK_SINGLE_CONTINUUM // Just track the total electron density. Computationally expensive to turn off. If off, tracks one electron distribution for each species (element) defined in input file, tracking the primary electrons that they release and the secondary electrons those electrons free from ALL species.
//#define NO_PLASMA // Disables all of the above. (i.e. primary - Photo,Auger,Fluoro - ionisation only)

/// Disable features
// Recommend to leave off plotting for now, doesn't seem to work anymore, and isn't so important since simulation times aren't too large anymore (but it's such a nice feature...).
#define NO_PLOTTING // Turns off live saves of the free-electron energy distribution to _live_plot.png. Disables use of python. 
//#define NO_BACKUP_SAVING // Disables the hourly saves of the data to  
#define TRACK_SINGLE_CONTINUUM // Just track the total electron density. Computationally expensive to turn off. If off, tracks one electron distribution for each species (element) defined in input file, tracking the primary electrons that they release and the secondary electrons those electrons free from ALL species.

/// Asynchronous solver
//#define NO_MINISTEPS   // Disables the asynchronous implementation of the solver, stepping the free (E-E) and bound (everything else) solvers together.
#define NO_MINISTEP_UPDATING  // ministep updates not working atm, possibly as lagrange polynomial is a bad thing to use here and I should just use linear interpolation.


// Dynamic grid modifiers
#define FIND_INITIAL_DIRAC // At point of first grid update (TODO make earlier) restart simulation, but using Dirac regions determined from the grid update. 

/// Dynamic grid debugging
//#define DEBUG_BOUND  // Turns on testing for negative bound state probabilities. Should narrow down cause (unless one is facing a dreaded memory corruption bug). Probably should be removed as this is no longer an issue.
//#define SWITCH_OFF_ALL_DYNAMIC_UPDATES  // The initial grid (spline knot basis) is used for the whole simulation. 
//#define SWITCH_OFF_DYNAMIC_BOUNDS // Disables updates of the MB energy bounds, the photoelectron energy bounds, and the transition energy. Should be equivalent to SWITCH_OFF_ALL_DYNAMIC_UPDATES if everything is working  


/// Analytical testing
//#define NO_ELECTRON_SOURCE  // Note this is equivalent to removing any #ELECTRON_SOURCE parameters in the .mol file.


/// Thread capping. Often an easier option than regenerating batch files. (If multiple are left on, uses the lowest cap)
#define THREAD_MAX_28
//#define THREAD_MAX_24
//#define THREAD_MAX_20
//#define THREAD_MAX_16
//#define THREAD_MAX_12
//#define THREAD_MAX_8
//#define THREAD_MAX_4
//#define THREAD_MAX_1

/// Experimental options for developers.

// #define NO_MIN_DIRAC_DENSITY // If this is off, peaks have to be a certain factor above the algorithmically determined transition point, so long as that point is below 600 eV. Good for binomial free electron distribution.  (Note to developers: currently that implementation could be improved (perhaps allowing this macro to be avoided) if transition energy is separated from cutoff energy, and made to be between the light auger and photopeaks.)


//#define BOUND_GD_HACK // A very primitive (unfinished) implementation of a bound transport correction. Not confirmed to be working as intended.

//#define CLASSIC_MANUAL_GRID // Implements an outdated static grid implementation for purposes of comparison.

//#define INFINITE_COULOMBLOG_CUTOFF  //Not confirmed to work with dynamic grid. Set cutoff for calculation of temperature/coulomb log to be infinite.

/////////////
#ifdef NO_MINISTEPS 
    #ifndef NO_MINISTEP_UPDATING
    #define NO_MINISTEP_UPDATING
    #endif
#endif //NO_MINISTEPS

#ifdef NO_PLASMA
    #ifndef NO_TBR
    #define NO_TBR
    #endif
    #ifndef NO_EE
    #define NO_EE
    #endif
    #ifndef NO_EII
    #define NO_EII
    #endif    
#endif // NO_PLASMA
#ifndef NO_PLOTTING // There's a bug with this right now
    #define NO_PLOTTING 
#endif
