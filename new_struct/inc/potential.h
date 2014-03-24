#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "coordsys.h"
#include "GSLInterface.h"

class Potential_JS{
	/* General base class for Potential_JS */
	public:
		inline virtual double Phi(VecDoub x){
			/* potential at Cartesian x 				 */
			/* -- to be overridden by derived classes -- */
			return 0.;
		}
		inline virtual VecDoub Forces(VecDoub x){
			/* forces at Cartesian x 					 */
			/* -- to be overridden by derived classes -- */
			return VecDoub(0.);
		}
		inline virtual double Vc(double R){
			/* circular velocity  at polar R 			*/
			/* -- to be overridden by derived classes -- */
			return 0.;
		}
		inline double H(VecDoub xv){
			/* returns Hamiltonian at Cartesian (x,v) */
			VecDoub x = {xv[0],xv[1],xv[2]};
			return 0.5*(xv[3]*xv[3]+xv[4]*xv[4]+xv[5]*xv[5])+Phi(x);
		}
		inline double Lz(VecDoub xv){
			/* returns z-component of angular momentum at Cartesian (x,v) */
			return xv[0]*xv[4]-xv[1]*xv[3];
		}
};

class StackelProlate_PerfectEllipsoid: public Potential_JS{
	/* Axisymmetric Stackel potential with perfect ellipsoidal density */
	/* Requires an oblate spheroidal coordinate system instance		   */
	private:
		double rho0, Const;
		OblateSpheroidCoordSys *CS;
	public:
		inline double alpha(){return CS->alpha();}
		inline double gamma(){return CS->gamma();}
		double G(double tau);
		double GPrime(double tau);
		VecDoub Vderivs(VecDoub tau);
		StackelProlate_PerfectEllipsoid(double Rho0, double alpha){
			CS = new OblateSpheroidCoordSys(alpha);
			Const = 2.*PI*Grav*Rho0*(-CS->alpha());
		}
		virtual ~StackelProlate_PerfectEllipsoid(){delete CS;}
		inline VecDoub x2tau(VecDoub x){return CS->x2tau(x);}
		inline VecDoub xv2tau(VecDoub x){return CS->xv2tau(x);}
		double Phi(VecDoub x);
		VecDoub Forces(VecDoub x);
		double Phi_tau(VecDoub tau);
		VecDoub x2ints(VecDoub x, VecDoub *tau = NULL);
};

class StackelTriaxial: public Potential_JS{
	/* Triaxial Stackel potential with perfect ellipsoidal density */
	/* Requires a confocal ellipsoidal coordinate system instance  */
	private:
		double rho0, Const, a, b, c, l, Flm, Elm, sinm;
		ConfocalEllipsoidalCoordSys *CS;	
		VecDoub Vderivs(VecDoub tau);
	public:
		StackelTriaxial(double Rho0, double alpha, double beta){
			CS = new ConfocalEllipsoidalCoordSys(alpha,beta);
			a = sqrt(-CS->alpha());b = sqrt(-CS->beta());c = sqrt(-CS->gamma());
			l = acos(c/a);double sinL = (1-(c/a)*(c/a));
			sinm = sqrt((1.-(b/a)*(b/a))/sinL);
			Flm = ellint_first(l,sinm);
			Elm = ellint_second(l,sinm);
			Const=2.*PI*Grav*Rho0*b*c/sqrt(sinL);
			// std::cout<<(-Const/a/a/sinL/sinm/sinm/3.*((1-sinm*sinm)*Flm+(2.*sinm*sinm-1.)*Elm-b*sinm*sinm/a*sqrt(sinL)*cos(l)))
			// 	<<" "<<-Const/a/a/sinL/sinm/sinm/(1-sinm*sinm)*(Elm-(1-sinm*sinm)*Flm-c*sinm*sinm*sin(l)/b)
			// 	<<" "<<-Const*0.5/a/a/3./sinL/(1-sinm*sinm)*(-(1-sinm*sinm)*Flm+(2.*(1-sinm*sinm)-1)*Elm+b*sin(l)*(b*b/c/c-(1-sinm*sinm))/c)<<std::endl;
		}
		virtual ~StackelTriaxial(){delete CS;}
		inline double alpha(){return CS->alpha();}
		inline double beta(){return CS->beta();}
		inline double gamma(){return CS->gamma();}
		double G(double tau);double GPrime(double tau);
		inline VecDoub x2tau(VecDoub x){return CS->x2tau(x);}
		inline VecDoub xv2tau(VecDoub x){return CS->xv2tau(x);}
		double Phi(VecDoub x);
		VecDoub Forces(VecDoub x);
		double Phi_tau(VecDoub tau);
		VecDoub tau2ints(VecDoub tau);
};

class Logarithmic: public Potential_JS{
	/* Two parameter axisymmetric logarithmic potential 	*/
	/* Phi = Vc^2/2 log(R^2+z^2/q^2)						*/
	private:
		double Vc2, q2;
	public:
		Logarithmic(double VC, double Q): Vc2(VC*VC), q2(Q*Q){}
		double Phi(VecDoub x);
		VecDoub Forces(VecDoub x);
};

class MiyamotoNagai: public Potential_JS{
	/* Axisymmetric Miyamoto-Nagai potential 	*/
	/* Phi = -GM/sqrt(R^2+(A+sqrt(z^2+b^2))^2)	*/
	private:
		double A, GM, Bq;
	public:
		MiyamotoNagai(double a, double gm, double b): A(a), GM(gm), Bq(b*b){}
		double Phi(VecDoub x);
		double Vc(double R);
		VecDoub Forces(VecDoub x);
};

class NFW: public Potential_JS{
	/* Axisymmetric NFW potential 								*/
	/* Phi = -GM*log(1.+sqrt(R*R+z*z/q2)/rs)/sqrt(R*R+z*z/q2)	*/
	private:
		double rs, GM, q2;
	public:
		NFW(double RS, double gm, double q): rs(RS), GM(gm), q2(q*q){}
		double Phi(VecDoub x);
		VecDoub Forces(VecDoub x);
};

// class GalPot: public Potential_JS{
// 	/* Wrapper for potentials produced by the GalPot code */
// 	private:
// 		GalaxyPotential *PhiWD;
// 	public:
// 		GalPot(std::string TpotFile);
// 		double Phi(VecDoub x);
// 		VecDoub Forces(VecDoub x);
// 		VecDoub getfreqs(double R);
// };

#endif
