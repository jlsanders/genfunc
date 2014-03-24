#ifndef ORBIT_H
#define ORBIT_H

#include "GSLInterface/GSLInterface.h"

class Orbit{
	/* class to integrate a phase-space point in a given Potential_JS */
	private:
		ode *O;
		std::vector<std::string> labels;
	public:
		std::vector<VecDoub> results;
		Orbit(Potential_JS *P, double eps = 1e-7);
		~Orbit(){delete O;}
		VecDoub integrate(VecDoub x, double t_interval, double step_size);
		void plot(int i, int j, std::string name = "orbit");
};
#endif