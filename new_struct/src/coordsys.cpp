/*======================================*/
/* 			Coordinate Systems 			*/
/*======================================*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_poly.h>
#include "utils.h"
#include "coordsys.h"

// ==========================================================================================
// Oblate Spheroidal Coordinate System
// ==========================================================================================

VecDoub OblateSpheroidCoordSys::x2tau(VecDoub x){
	/* Calculates tau given Cartesian x */
 	double R2 = x[0]*x[0]+x[1]*x[1];
 	double b = Alpha+Gamma-R2-x[2]*x[2];
 	double c = Alpha*Gamma-Gamma*R2-Alpha*x[2]*x[2];
 	VecDoub tau(3);
 	tau[0] = 0.5*(-b+sqrt(b*b-4.0*c)); 	/* lambda		*/
 	tau[2] = 0.5*(-b-sqrt(b*b-4.0*c)); 	/* nu			*/
 	tau[1] = atan2(x[1],x[0]); 			/* in radians	*/
 	return tau;
}

VecDoub OblateSpheroidCoordSys::xv2tau(VecDoub x){
	/* Calculates tau & tau_dot given Cartesian (x,v) */
	double R2 = x[0]*x[0]+x[1]*x[1];
  	double b = Alpha+Gamma-R2-x[2]*x[2];
  	double c = Alpha*Gamma-Gamma*R2-Alpha*x[2]*x[2];
  	double bdot = -2.0*(x[0]*x[3]+x[1]*x[4]+x[2]*x[5]);
  	double cdot = -2.0*(Gamma*(x[0]*x[3]+x[1]*x[4])+Alpha*x[2]*x[5]);
  	VecDoub tau(6);
  	tau[0] = 0.5*(-b+sqrt(b*b-4.0*c)); 	/* lambda		*/
 	tau[2] = 0.5*(-b-sqrt(b*b-4.0*c)); 	/* nu			*/
 	tau[1] = atan2(x[1],x[0]); 			/* in radians	*/
	tau[3] = 0.5*(-bdot+(b*bdot-2.0*cdot)/sqrt(b*b-4.0*c));
	tau[5] = 0.5*(-bdot-(b*bdot-2.0*cdot)/sqrt(b*b-4.0*c));
	tau[4] = (x[0]*x[4]-x[3]*x[1])/R2;
	return tau;
}

VecDoub OblateSpheroidCoordSys::derivs(VecDoub x){
	/* Calculates tau & derivatives of tau wrt to R and z given Cartesian x */

	double Rsq=x[0]*x[0]+x[1]*x[1], zsq=x[2]*x[2], A=-Alpha-Gamma+Rsq+zsq, 
	B=Alpha*Gamma-Gamma*Rsq-Alpha*zsq, temp=sqrt(A*A-4.*B),
    lamda=0.5*(A+temp), nu=0.5*(A-temp), 
    dldR=sqrt(Rsq)*(1.+(Gamma-Alpha+Rsq+zsq)/temp),
    dldz=x[2]*(1.+(Alpha-Gamma+Rsq+zsq)/temp),
    dvdR= sqrt(Rsq)*(1.-(Gamma-Alpha+Rsq+zsq)/temp),
    dvdz= x[2]*(1.-(Alpha-Gamma+Rsq+zsq)/temp);

	VecDoub derivs = {lamda,nu,dldR, dldz, dvdR, dvdz};

	return derivs;
}
VecDoub OblateSpheroidCoordSys::tau2p(VecDoub tau){
	/* Calculates p(tau) given 6D tau vector */
	double P2 = (tau[0]-tau[2])/(4.*(tau[0]+Alpha)*(tau[0]+Gamma));
	double R2 = (tau[2]-tau[0])/(4.*(tau[2]+Alpha)*(tau[2]+Gamma));
	VecDoub p = {P2*tau[3],R2*tau[5]};
	return p;
}

// ==========================================================================================
// Confocal Ellipsoidal Coordinate System
// ==========================================================================================

VecDoub ConfocalEllipsoidalCoordSys::x2tau(VecDoub x){
	/* Calculates tau given Cartesian x */
	double A = (Alpha+Beta+Gamma-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
	double B = (Alpha*Beta+Alpha*Gamma+Beta*Gamma-(Beta+Gamma)*x[0]*x[0]-(Alpha+Gamma)*x[1]*x[1]-(Alpha+Beta)*x[2]*x[2]);
	double C = Alpha*Beta*Gamma-Beta*Gamma*x[0]*x[0]-Alpha*Gamma*x[1]*x[1]-Alpha*Beta*x[2]*x[2];
	double tau[3];		
	gsl_poly_solve_cubic (A,B,C,&tau[2],&tau[1],&tau[0]);
	VecDoub tau_v = {tau[0],tau[1],tau[2]};
	return tau_v;
}

VecDoub ConfocalEllipsoidalCoordSys::tau2x(VecDoub tau){
	/* Calculates Cartesian at tau */
	double x = sqrt((tau[0]+Alpha)*(tau[1]+Alpha)*(tau[2]+Alpha)/(Alpha-Gamma)/(Alpha-Beta));
	double y = sqrt((tau[0]+Beta)*(tau[1]+Beta)*(tau[2]+Beta)/(Beta-Gamma)/(Beta-Alpha));
	double z = sqrt((tau[0]+Gamma)*(tau[1]+Gamma)*(tau[2]+Gamma)/(Gamma-Beta)/(Gamma-Alpha));
	return {x,y,z};
}

VecDoub ConfocalEllipsoidalCoordSys::xv2tau(VecDoub x){
	/* Calculates tau & taudot given 6D Cartesian (x,v) */
	double A = (Alpha+Beta+Gamma-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
	double B = (Alpha*Beta+Alpha*Gamma+Beta*Gamma-(Beta+Gamma)*x[0]*x[0]-(Alpha+Gamma)*x[1]*x[1]-(Alpha+Beta)*x[2]*x[2]);
	double C = Alpha*Beta*Gamma-Beta*Gamma*x[0]*x[0]-Alpha*Gamma*x[1]*x[1]-Alpha*Beta*x[2]*x[2];
	double tau[3];

	gsl_poly_solve_cubic (A,B,C,&tau[2],&tau[1],&tau[0]);

	double Aup = x[0]*x[3]+x[1]*x[4]+x[2]*x[5];
	double Bup = (Beta+Gamma)*x[0]*x[3]+(Alpha+Gamma)*x[1]*x[4]+(Alpha+Beta)*x[2]*x[5];
	double Cup = Beta*Gamma*x[0]*x[3]+Alpha*Gamma*x[1]*x[4]+Alpha*Beta*x[2]*x[5];

	VecDoub tau_v(6); for(int i=0; i<3; i++) tau_v[i]=tau[i];
	tau_v[3] = ((2.0*tau[0]*tau[0]*Aup+2.0*tau[0]*Bup+2.0*Cup)/(3.0*tau[0]*tau[0]+2.0*tau[0]*A+B));
	tau_v[4] = ((2.0*tau[1]*tau[1]*Aup+2.0*tau[1]*Bup+2.0*Cup)/(3.0*tau[1]*tau[1]+2.0*tau[1]*A+B));
	tau_v[5] = ((2.0*tau[2]*tau[2]*Aup+2.0*tau[2]*Bup+2.0*Cup)/(3.0*tau[2]*tau[2]+2.0*tau[2]*A+B));
	return tau_v;
}
VecDoub ConfocalEllipsoidalCoordSys::tau2p(VecDoub tau){
	/* Calculates (P^2 l_dot^2,...)(tau) given 6D tau vector */
	double P2 = (tau[0]-tau[1])*(tau[0]-tau[2])/(4.*(tau[0]+Alpha)*(tau[0]+Beta)*(tau[0]+Gamma));
	double Q2 = (tau[1]-tau[2])*(tau[1]-tau[0])/(4.*(tau[1]+Alpha)*(tau[1]+Beta)*(tau[1]+Gamma));
	double R2 = (tau[2]-tau[0])*(tau[2]-tau[1])/(4.*(tau[2]+Alpha)*(tau[2]+Beta)*(tau[2]+Gamma));
	VecDoub pt = {P2*tau[3]*tau[3],Q2*tau[4]*tau[4],R2*tau[5]*tau[5]};
	return pt;
}

VecDoub ConfocalEllipsoidalCoordSys::derivs(VecDoub x){
	/* Finds the derivatives of tau coordinates wrt to Cartesian x */
	/* The vector returned has the tau coordinates as the first    */
	/* three elements and then dl/dx,dm/dx,dn/dx,dl/dy...		   */	
	VecDoub derivs = x2tau(x);

	VecDoub dA = {-2*x[0],-2*x[1],-2*x[2]};
	VecDoub dB = {-2*(Beta+Gamma)*x[0],-2*(Alpha+Gamma)*x[1],-2*(Alpha+Beta)*x[2]};
	VecDoub dC = {-2*Beta*Gamma*x[0],-2*Alpha*Gamma*x[1],-2*Alpha*Beta*x[2]};
	double A = (Alpha+Beta+Gamma-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
	double B = (Alpha*Beta+Alpha*Gamma+Beta*Gamma-(Beta+Gamma)*x[0]*x[0]-(Alpha+Gamma)*x[1]*x[1]-(Alpha+Beta)*x[2]*x[2]);
	VecDoub denom = {(3.*derivs[0]*derivs[0]+2.*A*derivs[0]+B),(3.*derivs[1]*derivs[1]+2.*A*derivs[1]+B),(3.*derivs[2]*derivs[2]+2.*A*derivs[2]+B)};

	// dtaudx
	derivs.push_back(-(dC[0]+derivs[0]*dB[0]+derivs[0]*derivs[0]*dA[0])/denom[0]);
	derivs.push_back(-(dC[0]+derivs[1]*dB[0]+derivs[1]*derivs[1]*dA[0])/denom[1]);
	derivs.push_back(-(dC[0]+derivs[2]*dB[0]+derivs[2]*derivs[2]*dA[0])/denom[2]);
	
	// dtaudy
	derivs.push_back(-(dC[1]+derivs[0]*dB[1]+derivs[0]*derivs[0]*dA[1])/denom[0]);
	derivs.push_back(-(dC[1]+derivs[1]*dB[1]+derivs[1]*derivs[1]*dA[1])/denom[1]);
	derivs.push_back(-(dC[1]+derivs[2]*dB[1]+derivs[2]*derivs[2]*dA[1])/denom[2]);

	// dtaudz
	derivs.push_back(-(dC[2]+derivs[0]*dB[2]+derivs[0]*derivs[0]*dA[2])/denom[0]);
	derivs.push_back(-(dC[2]+derivs[1]*dB[2]+derivs[1]*derivs[1]*dA[2])/denom[1]);
	derivs.push_back(-(dC[2]+derivs[2]*dB[2]+derivs[2]*derivs[2]*dA[2])/denom[2]);

	return derivs;
}

// ==========================================================================================
// coordsys.cpp