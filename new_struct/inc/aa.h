#ifndef AA_H
#define AA_H

#include "potential.h"

// ==========================================================================================
// Axisymmetric Stackel Perfect Ellipsoid Potential
// ==========================================================================================

struct root_struct_axi{
	StackelProlate_PerfectEllipsoid *P;
	VecDoub Ints;
	root_struct_axi(StackelProlate_PerfectEllipsoid *PP, VecDoub ints)
		:P(PP),Ints(ints){}
};

struct action_struct_axi{
	root_struct_axi RS;
	double taubargl, Deltagl;
	action_struct_axi(StackelProlate_PerfectEllipsoid *PP, VecDoub ints, double tb, double Dl)
		:RS(PP,ints),taubargl(tb),Deltagl(Dl){} 
};

class Actions_AxisymmetricStackel{
	private:
		StackelProlate_PerfectEllipsoid *Pot;
		VecDoub find_limits(VecDoub x, VecDoub ints);
	public:
		Actions_AxisymmetricStackel(StackelProlate_PerfectEllipsoid *pot): Pot(pot){};
		VecDoub actions(VecDoub x);
};

// ==========================================================================================
// Triaxial Stackel Perfect Ellipsoid Potential
// ==========================================================================================

struct root_struct_triax{
	StackelTriaxial *P;
	VecDoub Ints;
	root_struct_triax(StackelTriaxial *PP, VecDoub ints)
		:P(PP),Ints(ints){}
};

struct action_struct_triax{
	StackelTriaxial *P;
	VecDoub Ints;
	double taubargl, Deltagl;
	action_struct_triax(StackelTriaxial *PP, VecDoub ints, double tb, double Dl)
		:P(PP),Ints(ints),taubargl(tb),Deltagl(Dl){} 
};

class Actions_TriaxialStackel{
	private:
		StackelTriaxial *Pot;
		VecDoub find_limits(VecDoub x,VecDoub ints);
	public:
		Actions_TriaxialStackel(StackelTriaxial *pot): Pot(pot){};
		VecDoub actions(VecDoub x, bool freq_yes=0);	
};

#endif
// ==========================================================================================