/*======================================*/
/* 			    Potential_JS 			*/
/*======================================*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "coordsys.h"
#include "coordtransforms.h"
#include "potential.h"

// ==========================================================================================
// Prolate Stackel Perfect Ellipsoid Potential
// ==========================================================================================

double StackelProlate_PerfectEllipsoid::G(double tau){
	/* de Zeeuw's G(tau) function */
	double Gamma = CS->gamma(), sqG =sqrt(-Gamma/(tau+Gamma));
	return Const*sqG*atan(1./sqG);
}
double StackelProlate_PerfectEllipsoid::GPrime(double tau){
	/* derivative of G wrt tau */
	double Gamma = CS->gamma(), sqG =sqrt(-Gamma/(tau+Gamma));
	return 0.5*Const*sqG*sqG*(sqG*atan(1./sqG)/Gamma+1./tau);
}

VecDoub StackelProlate_PerfectEllipsoid::Vderivs(VecDoub tau){
	/* Calculates the derivatives of the Potential_JS wrt tau given tau=(l,v) */
	double Gl = G(tau[0]), Gv = G(tau[1]), Gamma = CS->gamma();
	double dVdl = 	(-GPrime(tau[0])*(tau[0]+Gamma)-Gl
			+(Gl*(tau[0]+Gamma)-Gv*(tau[1]+Gamma))/(tau[0]-tau[1]))
			/(tau[0]-tau[1]);
	double dVdv = 	(GPrime(tau[1])*(tau[1]+Gamma)+Gv
			-(Gl*(tau[0]+Gamma)-Gv*(tau[1]+Gamma))/(tau[0]-tau[1]))
			/(tau[0]-tau[1]);
	VecDoub derivs = {dVdl, dVdv};
	return derivs;
}


VecDoub StackelProlate_PerfectEllipsoid::Forces(VecDoub x){
	/* forces at Cartesian x */
	VecDoub derivs = CS->derivs(x);
	VecDoub tau = {derivs[0],derivs[1]};
	VecDoub Vderiv = Vderivs(tau);

	double dvdR = -Vderiv[0]*derivs[2]-Vderiv[1]*derivs[4];
	double R = sqrt(x[0]*x[0]+x[1]*x[1]);
	VecDoub result ={ 	x[0]*dvdR/R, x[1]*dvdR/R,
			  			-Vderiv[0]*derivs[3]-Vderiv[1]*derivs[5]};
	return result;
}

double StackelProlate_PerfectEllipsoid::Phi_tau(VecDoub tau){
	/* Potential at tau */
	double Gamma = CS->gamma();
	return -((tau[0]+Gamma)*G(tau[0])-(tau[2]+Gamma)*G(tau[2]))/(tau[0]-tau[2]);
}

double StackelProlate_PerfectEllipsoid::Phi(VecDoub x){
	/* potential at Cartesian x */
	VecDoub tau = CS->x2tau(x);
	return Phi_tau(tau);
}

VecDoub StackelProlate_PerfectEllipsoid::x2ints(VecDoub x, VecDoub *tau){
	VecDoub Ints = {H(x), 0.5*pow(Lz(x),2.)};
	if(!tau) (*tau) = CS->xv2tau(x);
	Ints.push_back(((*tau)[0]+CS->gamma())*(Ints[0]-(Ints[1]/((*tau)[0]+CS->alpha()))
		+G((*tau)[0]))-(pow(((*tau)[3]*((*tau)[0]-(*tau)[2])),2.0))/(8.0*((*tau)[0]+CS->alpha())
		*((*tau)[0]+CS->gamma())));
	printVector(Ints);
	//Ints.push_back(tau[0]); Ints.push_back(tau[2]);
	return Ints;
}

// ==========================================================================================
// Triaxial Stackel Perfect Ellipsoid Potential
// ==========================================================================================

double taugl, acgl, bcgl, c2gl;
double G_integrand(double s,void*){
	return sqrt(c2gl+bcgl*s*s)/sqrt(c2gl+acgl*s*s)/(c2gl+(taugl-c2gl)*s*s);
}

double StackelTriaxial::G(double tau){
	/* de Zeeuw's perfect ellipsoid G function 			*/
	/* calculated using GL integration of eq B9 of dZ85 */
	taugl=tau;c2gl=c*c;acgl=a*a-c2gl; bcgl=b*b-c2gl;
	return Const*sin(l)*a*c*GaussLegendreQuad(&G_integrand,0.,1.);
}

double GP_integrand(double s, void*){
	double p = c2gl+(taugl-c2gl)*s*s;
	return -sqrt(c2gl+bcgl*s*s)/sqrt(c2gl+acgl*s*s)*s*s/p/p;
}
double StackelTriaxial::GPrime(double tau){
	/* de Zeeuw's perfect ellipsoid GPrime function 	*/
	/* calculated using GL integration of eq B9 of dZ85 */
	taugl=tau;c2gl=c*c;acgl=a*a-c2gl;bcgl=b*b-c2gl;
	return Const*sin(l)*a*c*GaussLegendreQuad(&GP_integrand,0.,1.);
}

/*
// Attempts at analytic expressions for above -- doesn't work

double StackelTriaxial::G(double tau){
	// std::cout<<l<<" "<<sinm<<" "<<-(tau+CS->alpha())/(CS->alpha()-CS->gamma())<<std::endl;
	// std::cout<<ellint_third(l,sinm,-(tau+CS->alpha())/(CS->alpha()-CS->gamma()))<<std::endl;
	std::cout<<ellint_third(l,sinm,0.)<<" "<<Flm<<" "<<Elm<<std::endl;
	return Const/(tau+CS->alpha())*
		   ((tau+CS->beta())*ellint_third(l,sinm,-(tau+CS->alpha())/(CS->alpha()-CS->gamma()))
		   +(CS->alpha()-CS->beta())*Flm);
}

double StackelTriaxial::GPrime(double tau){
	return 0.5*Const/(tau+CS->alpha())/(tau+CS->gamma())*
		   (((tau+CS->alpha())*(tau+CS->gamma())-(tau+CS->alpha())*(tau+CS->beta())-(tau+CS->gamma())*(tau+CS->beta()))
		   	*ellint_third(l,sinm,-(tau+CS->alpha())/(CS->alpha()-CS->gamma()))/(tau+CS->alpha())
		   +((tau+CS->alpha())*(CS->beta()-CS->gamma())+(tau+CS->gamma())*(CS->beta()-CS->alpha()))
		    *Flm/(tau+CS->alpha())
		   -(tau+CS->alpha())*b*c*sin(l)/CS->alpha()+(CS->gamma()-CS->alpha())*Elm);
}
*/

VecDoub StackelTriaxial::Vderivs(VecDoub tau){
	/* Calculates the derivatives of the Potential_JS wrt tau given tau=(l,m,v) */
	double Gl = G(tau[0]), Gm = G(tau[1]), Gn = G(tau[2]), Alpha = CS->alpha(), Gamma = CS->gamma();
	double  F_l = (tau[0]+Alpha)*(tau[0]+Gamma)*Gl,
			F_m = (tau[1]+Alpha)*(tau[1]+Gamma)*Gm,
			F_n = (tau[2]+Alpha)*(tau[2]+Gamma)*Gn;
	double dlm = tau[0]-tau[1], dln = tau[0]-tau[2];
	double dml = tau[1]-tau[0], dmn = tau[1]-tau[2];
	double dnl = tau[2]-tau[0], dnm = tau[2]-tau[1];

	double dVdl = -F_l*GPrime(tau[0])/Gl/dlm/dln-(2.*tau[0]+Alpha+Gamma)*Gl/dlm/dln+F_l/dlm/dln*(1./dln+1./dlm)
					-F_m/dmn/dml/dml-F_n/dnl/dnl/dnm;
	double dVdm = -F_m*GPrime(tau[1])/Gm/dml/dmn-(2.*tau[1]+Alpha+Gamma)*Gm/dmn/dml+F_m/dml/dmn*(1./dmn+1./dml)
					-F_l/dlm/dlm/dln-F_n/dnl/dnm/dnm;
	double dVdn = -F_n*GPrime(tau[2])/Gn/dnl/dnm-(2.*tau[2]+Alpha+Gamma)*Gn/dnm/dnl+F_n/dnl/dnm*(1./dnm+1./dnl)
					-F_l/dlm/dln/dln-F_m/dmn/dmn/dml;
	VecDoub derivs = {dVdl, dVdm, dVdn};
	return derivs;
}

VecDoub StackelTriaxial::Forces(VecDoub x){
	/* forces at Cartesian x */
	VecDoub derivs = CS->derivs(x);
	VecDoub tau = {derivs[0],derivs[1],derivs[2]};
	VecDoub Vderiv = Vderivs(tau);
	VecDoub result ={ 	-Vderiv[0]*derivs[3]-Vderiv[1]*derivs[4]-Vderiv[2]*derivs[5],
						-Vderiv[0]*derivs[6]-Vderiv[1]*derivs[7]-Vderiv[2]*derivs[8],
			  			-Vderiv[0]*derivs[9]-Vderiv[1]*derivs[10]-Vderiv[2]*derivs[11]};
	return result;
}

double StackelTriaxial::Phi_tau(VecDoub tau){
	/* Potential at tau */
	double  F_l = (tau[0]+CS->alpha())*(tau[0]+CS->gamma())*G(tau[0]),
			F_m = (tau[1]+CS->alpha())*(tau[1]+CS->gamma())*G(tau[1]),
			F_n = (tau[2]+CS->alpha())*(tau[2]+CS->gamma())*G(tau[2]);
	return -F_l/(tau[0]-tau[1])/(tau[0]-tau[2])-F_m/(tau[1]-tau[2])/(tau[1]-tau[0])-F_n/(tau[2]-tau[0])/(tau[2]-tau[1]);
}

double StackelTriaxial::Phi(VecDoub x){
	/* potential at Cartesian x */
	VecDoub tau = CS->x2tau(x);
	return Phi_tau(tau);
}

VecDoub StackelTriaxial::tau2ints(VecDoub tau){
	// Needs correct equations
	VecDoub pp = CS->tau2p(tau);
	double X = 0.5*pp[0]-(tau[0]+CS->alpha())*(tau[0]+CS->gamma())*G(tau[0])/(tau[0]-tau[1])/(tau[0]-tau[2]);
	double Y = 0.5*pp[1]-(tau[1]+CS->alpha())*(tau[1]+CS->gamma())*G(tau[1])/(tau[1]-tau[0])/(tau[1]-tau[2]);
	double Z = 0.5*pp[2]-(tau[2]+CS->alpha())*(tau[2]+CS->gamma())*G(tau[2])/(tau[2]-tau[1])/(tau[2]-tau[0]);
	VecDoub Ints = {X+Y+Z};
	double J =(tau[1]+tau[2])*X+(tau[2]+tau[0])*Y+(tau[0]+tau[1])*Z;
	double K = tau[1]*tau[2]*X+tau[2]*tau[0]*Y+tau[0]*tau[1]*Z;
	Ints.push_back((CS->alpha()*(CS->alpha()*Ints[0]+J)+K)/(CS->alpha()-CS->gamma()));	
	Ints.push_back((CS->gamma()*(CS->gamma()*Ints[0]+J)+K)/(CS->gamma()-CS->alpha()));
	return Ints;
}

// ==========================================================================================
// Logarithmic Potential
// ==========================================================================================

double Logarithmic::Phi(VecDoub x){
	/* potential at Cartesian x */
	return Vc2/2.*log(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]/q2);
}
VecDoub Logarithmic::Forces(VecDoub x){
	/* Forces at Cartesian x */
	double r = Vc2/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]/q2);
	VecDoub F = {-x[0]*r,-x[1]*r,-x[2]*r/q2};
	return F;
}

// ==========================================================================================
// Miyamoto-Nagai Potential
// ==========================================================================================

double MiyamotoNagai::Phi(VecDoub x){
	/* potential at Cartesian x */
	double AZB=A+sqrt(x[2]*x[2]+Bq);
	return -GM/sqrt(x[0]*x[0]+x[1]*x[1]+AZB*AZB);
}

double MiyamotoNagai::Vc(double R){
	double t = A+sqrt(Bq);
	double F = 1./(R*R+t*t);
	return sqrt(GM*R*R*F*sqrt(F));
}

VecDoub MiyamotoNagai::Forces(VecDoub x){
	/* Forces at Cartesian x */
	double  ZB=sqrt(x[2]*x[2]+Bq), AZB = A + ZB;
	double f = 1./(x[0]*x[0]+x[1]*x[1]+AZB*AZB),rtF=sqrt(f);
	VecDoub Force = {-GM*f*rtF*x[0], -GM*x[1]*f*rtF, ZB? -GM*x[2]*f*rtF*AZB/ZB : 0.};
	return Force;
}

// ==========================================================================================
// NFW Potential
// ==========================================================================================

double NFW::Phi(VecDoub x){
	/* potential at Cartesian x */
	double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]/q2);
	return -GM*log(1.+r/rs)/r;
}
VecDoub NFW::Forces(VecDoub x){
	/* Forces at Cartesian x */
	double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]/q2);
	double dpdr = GM*(log(1.+r/rs)/r-1./rs/(1.+r/rs))/r;
	VecDoub Force = {-dpdr*x[0]/r, -dpdr*x[1]/r, -dpdr*x[2]/r/q2};
	return Force;
}

// ==========================================================================================
// GalPot Potential for interface with Walter Dehnen's code
// ==========================================================================================

// GalPot::GalPot(std::string TpotFile){
// 	std::ifstream file;
// 	file.open(TpotFile);
// 	PhiWD=new GalaxyPotential(file);
// 	file.close();
// }

// double GalPot::Phi(VecDoub x){
// 	/* potential at Cartesian x */
// 	double R = sqrt(x[0]*x[0]+x[1]*x[1]);
// 	return conv::kpcMyr2kmsSq*(*PhiWD)(R,x[2]);
// }
// VecDoub GalPot::Forces(VecDoub x){
// 	/* Forces at Cartesian x */
// 	double R = sqrt(x[0]*x[0]+x[1]*x[1]);
// 	double dR, dz;
// 	(*PhiWD)(R,x[2],dR,dz);
// 	VecDoub f = {-conv::kpcMyr2kmsSq*x[0]*dR/R,-conv::kpcMyr2kmsSq*x[1]*dR/R, -conv::kpcMyr2kmsSq*dz};
// 	return f;
// }

// VecDoub GalPot::getfreqs(double R){
// 	/* frequencies: kappa, omega_c and nu at polar R */
// 	Frequencies freqs=(*PhiWD).KapNuOm(R);
// 	VecDoub freqs_v = {conv::kpcMyr2kms*freqs(0),conv::kpcMyr2kms*freqs(1),conv::kpcMyr2kms*freqs(2)};
// 	return freqs_v;
// }

// ==========================================================================================
// potential.cpp