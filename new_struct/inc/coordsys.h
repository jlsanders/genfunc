#ifndef COORDSYS_H
#define COORDSYS_H

#include "utils.h"

class OblateSpheroidCoordSys{
	private:
		double Alpha;
		const double Gamma;
	public:
		OblateSpheroidCoordSys(double alpha): Alpha(alpha), Gamma(-1.){}
		inline double alpha(void){ return Alpha;}
		inline double gamma(void){ return Gamma;}
		inline void newalpha(double a){ Alpha = a;}
		VecDoub x2tau(VecDoub x);
		VecDoub xv2tau(VecDoub x);
		VecDoub derivs(VecDoub x);
		VecDoub tau2p(VecDoub tau);
};


class ConfocalEllipsoidalCoordSys{
	private:
		double Alpha, Beta, Gamma;
	public:
		ConfocalEllipsoidalCoordSys(double a, double b): Alpha(a), Beta(b), Gamma(-1.){}
		inline double alpha(void){ return Alpha;}
		inline double beta(void){ return Beta;}
		inline double gamma(void){ return Gamma;}
		inline void newalpha(double a){ Alpha = a;}
		inline void newbeta(double b){ Beta = b;}
		VecDoub x2tau(VecDoub x);
		VecDoub tau2x(VecDoub tau);
		VecDoub xv2tau(VecDoub x);
		VecDoub tau2p(VecDoub tau);
		VecDoub derivs(VecDoub x);
};

#endif