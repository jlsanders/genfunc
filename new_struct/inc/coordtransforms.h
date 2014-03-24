#ifndef COORDTRANSFORMS_H
#define COORDTRANSFORMS_H

#include "utils.h"

namespace conv{

// Conventions
const VecDoub StandardSolar {8.0,0.0,11.1,232.24,7.25}; // in km/s with v_c = 220 km/s
const VecDoub StandardSolar2 {8.29,0.0,0.01135,0.257055,0.0074148}; // in kpc/Myr with v_c=239.1km/s
const VecDoub StandardSolarPAUL {8.29,0.0,11.1,251.34,7.25}; // in km/s with v_c = 239.1 km/s
const double masyr2radMyr = 4.8481368111e-3;
const double deg2rad = 0.017453292;
const double kpcMyr2kms = 977.775;
const double kpcMyr2kmsSq = kpcMyr2kms*kpcMyr2kms;
const double kms2kpcMyr = 1./kpcMyr2kms;
const double kms2kpcGyr = 1000.*kms2kpcMyr;
const double RA_GP=3.36603292,decGP=0.473477282,lCP=2.145566725;
const double PM_Const = 4.74057170372;

// ======================================================================================
// Cartesian <==> Polar 
VecDoub CartesianToPolar(const VecDoub& Cartesian);					
VecDoub PolarToCartesian(const VecDoub& Polar);

// ======================================================================================
// Galactic <==> Cartesian
VecDoub GalacticToCartesian(const VecDoub &Galactic,
								      const VecDoub& SolarPosition);

VecDoub GalacticToCartesian(const VecDoub &Galactic);

VecDoub CartesianToGalactic(const VecDoub &Cartesian,
									const VecDoub& SolarPosition);

VecDoub CartesianToGalactic(const VecDoub &Cartesian);
// ======================================================================================
// Galactic <==> Polar 
VecDoub PolarToGalactic(const VecDoub &Polar,
									const VecDoub& SolarPosition);
VecDoub PolarToGalactic(const VecDoub &Polar);

VecDoub GalacticToPolar(const VecDoub &Galactic,
									const VecDoub& SolarPosition);
	
VecDoub GalacticToPolar(const VecDoub &Galactic);
// ======================================================================================
// Equatorial <==> Galactic 
VecDoub EquatorialToGalactic(const VecDoub &Eq);
std::vector<VecDoub> EquatorialToGalacticwithErrors(const VecDoub &E,const VecDoub &EqE);
VecDoub GalacticToEquatorial(const VecDoub &Galactic);
}
#endif
