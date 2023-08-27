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
//#define NO_EE   // Electron-electron scattering. This can break the dynamic grid late in the simulation at low energies.
//#define NO_EII    // Electron impact ionisation
//--

//#define NO_BACKUP_SAVING // Disables the hourly saves of the data to output/backup_data/


// Asynchronous solver
#define NO_MINISTEP_UPDATING  // ministep updates not working atm, possibly as lagrange polynomial is a bad thing to use here and I should just use linear interpolation.


// Dynamic grid debugging
//#define DEBUG_BOUND  // Turns on testing for negative bound state probabilities. Should narrow down cause (unless one is facing a dreaded memory corruption bug). Probably should be removed as this is no longer an issue.
//#define SWITCH_OFF_ALL_DYNAMIC_UPDATES  // The initial grid (spline knot basis) is used for the whole simulation. 
//#define SWITCH_OFF_DYNAMIC_BOUNDS // Disables updates of the MB energy bounds, the photoelectron energy bounds, and the transition energy. Should be equivalent to SWITCH_OFF_ALL_DYNAMIC_UPDATES if everything is working  




/// Experimental options for developers.


//#define BOUND_GD_HACK // A very primitive (unfinished) implementation of a bound transport correction. Not confirmed to be working as intended.

//#define CLASSIC_MANUAL_GRID // Implements an outdated static grid implementation for purposes of comparison.

//#define INFINITE_COULOMBLOG_CUTOFF  //Not confirmed to work with dynamic grid. Set cutoff for calculation of temperature/coulomb log to be infinite.
/////////////


#ifdef NO_MINISTEPS 
    #ifndef NO_MINISTEP_UPDATING
    #define NO_MINISTEP_UPDATING
    #endif
#endif
