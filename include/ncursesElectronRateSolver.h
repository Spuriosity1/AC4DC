#include "ElectronRateSolver.h"
#ifdef PYBIND
#include "Plotting.h"
#endif //PYBIND
#include "Display.h"
#include <fstream>

class ncursesElectronRateSolver : public ElectronRateSolver {
	public:
		ncursesElectronRateSolver(const char* filename, ofstream& log) :
			ElectronRateSolver(filename, log){};

	void pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update);
	int post_ode_step(ofstream& _log, size_t& n);
	protected:

	#ifdef PYBIND
    // Display stuff
    Plotting py_plotter;
	#endif //PYBIND
};
