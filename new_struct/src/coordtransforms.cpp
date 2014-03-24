/*==============================================*/
/* 			Coordinate Transformations 			*/
/*==============================================*/

// ======================================================================================
#include <vector>
#include <cmath>
#include <iostream>
#include "utils.h"
#include "coordtransforms.h"

namespace conv{
// ======================================================================================
// Cartesian <==> Polar 
VecDoub CartesianToPolar(const VecDoub& Cartesian){
	// X,Y,Z -> R,phi,z
	VecDoub Polar {	sqrt(Cartesian[0]*Cartesian[0]+Cartesian[1]*Cartesian[1]),
					atan2(Cartesian[1],Cartesian[0]),
					Cartesian[2]};
	if(Cartesian.size()==3)	return Polar;
	// vx,vy,vz -> vR,vphi,vz
	else{
		double cp = cos(Polar[1]), sp = sin(Polar[1]);
		VecDoub PolarVel {	Cartesian[3]*cp+Cartesian[4]*sp,Cartesian[4]*cp-Cartesian[3]*sp,
					        Cartesian[5]};
		for (	VecDoub::iterator it = PolarVel.begin();
				it != PolarVel.end(); ++it) Polar.push_back(*it);
		return Polar;
		}
}
							
VecDoub PolarToCartesian(const VecDoub& Polar){
	// R,phi,z -> X,Y,Z
	double cp = cos(Polar[1]), sp = sin(Polar[1]);
	VecDoub Cartesian {	Polar[0]*cp,
						Polar[0]*sp,
						Polar[2]};
	if(Polar.size()==3) return Cartesian;
	// vR,vphi,vz -> vx,vy,vz
	else{
		VecDoub CartVel {Polar[3]*cp-Polar[4]*sp,Polar[4]*cp+Polar[3]*sp,Polar[5]};
		for (	VecDoub::iterator it = CartVel.begin();
				it != CartVel.end(); ++it) Cartesian.push_back(*it);
		return Cartesian;
		}
}

// ======================================================================================
// Galactic <==> Cartesian
VecDoub GalacticToCartesian(const VecDoub &Galactic,
								      const VecDoub& SolarPosition){
	// l,b,s->X,Y,Z
	double 	cl = cos(Galactic[0]), sl = sin(Galactic[0]),
			cb = cos(Galactic[1]), sb = sin(Galactic[1]);
	VecDoub Cartesian {	SolarPosition[0]-Galactic[2]*cb*cl,
						-Galactic[2]*cb*sl,
						Galactic[2]*sb+SolarPosition[1]};

	if(Galactic.size()==3)return Cartesian;
	// vlos,mu_lcos(b),mu_b -> vx,vy,vz
	// in units km/s, rad/Myr -> km/s
	else{ 	double vl = PM_Const*Galactic[2]*Galactic[4];double vb = PM_Const*Galactic[2]*Galactic[5];
			double tmp = cb*Galactic[3]-sb*vb;
			double vx = cl*tmp-sl*vl+SolarPosition[2];
			double vy = sl*tmp+cl*vl+SolarPosition[3];
			double vz = sb*Galactic[3]+cb*vb+SolarPosition[4];
			VecDoub CartVel{-vx,-vy,vz};
	  		for (	VecDoub::iterator it = CartVel.begin();
					it != CartVel.end(); ++it) Cartesian.push_back(*it);
				return Cartesian;
	}
}

VecDoub GalacticToCartesian(const VecDoub &Galactic){
  return GalacticToCartesian(Galactic,StandardSolar);
}

VecDoub CartesianToGalactic(const VecDoub &Cartesian,
									const VecDoub& SolarPosition){
	// X,Y,Z->l,b,s
	double tmp1 = SolarPosition[0]-Cartesian[0];
	double tmp2 = -Cartesian[1]; double tmp3 = Cartesian[2]-SolarPosition[1];
	double Distance = sqrt(tmp1*tmp1+tmp2*tmp2+tmp3*tmp3);
	VecDoub Galactic {	atan2(tmp2,tmp1),
						asin(tmp3/Distance),
						Distance};
	if(Cartesian.size()==3)return Galactic;
	// vx,vy,vz -> vlos,mu_lcos(b),mu_b
	// in units km/s -> km/s rad/Myr
	else{ 	double vx=-Cartesian[3]-SolarPosition[2];
			double vy = -Cartesian[4]-SolarPosition[3];
			double vz = Cartesian[5]-SolarPosition[4];
			double 	cl = cos(Galactic[0]), sl = sin(Galactic[0]),
			cb = cos(Galactic[1]), sb = sin(Galactic[1]);
			VecDoub GalVel {vx*cl*cb+vy*sl*cb+vz*sb,(-vx*sl+vy*cl)/(PM_Const*Distance),
				        	(-vx*cl*sb-vy*sl*sb+vz*cb)/(PM_Const*Distance)};
			for (	VecDoub::iterator it = GalVel.begin();
					it != GalVel.end(); ++it) Galactic.push_back(*it);
			return Galactic;
		}
}

VecDoub CartesianToGalactic(const VecDoub &Cartesian){
  return CartesianToGalactic(Cartesian,StandardSolar);
}
// ======================================================================================
// Galactic <==> Polar 
VecDoub PolarToGalactic(const VecDoub &Polar,
									const VecDoub& SolarPosition){
	// R,phi,z -> l,b,s; SolarPosition=(R0,z0)
	VecDoub Cart = PolarToCartesian(Polar);
	VecDoub Galactic = CartesianToGalactic(Cart,SolarPosition);
	return Galactic;
	}
	
VecDoub PolarToGalactic(const VecDoub &Polar){
  return PolarToGalactic(Polar,StandardSolar);
}

VecDoub GalacticToPolar(const VecDoub &Galactic,
									const VecDoub& SolarPosition){
	// l,b,s->R,phi,z; SolarPosition=(R0,z0)
	VecDoub Cart = GalacticToCartesian(Galactic,SolarPosition);
	VecDoub Polar = CartesianToPolar(Cart);
	return Polar;
	}
	
VecDoub GalacticToPolar(const VecDoub &Galactic){
  return GalacticToPolar(Galactic,StandardSolar);
}
// ======================================================================================
// Equatorial <==> Galactic 

VecDoub EquatorialToGalactic(const VecDoub &Equatorial){
	//alpha, dec, s => l,b,s
	double alpha = Equatorial[0], delta = Equatorial[1];
	double cd = cos(delta), sd = sin(delta);
	double b=asin(sin(decGP)*sd+cos(decGP)*cd*cos(alpha-RA_GP));
	double l=lCP-atan2(cd*sin(alpha-RA_GP),cos(decGP)*sd-sin(decGP)*cd*cos(alpha-RA_GP));
	if(l<0.)l+=2.*PI;
	VecDoub Galactic {l,b,Equatorial[2]};
	if(Equatorial.size()==3)return Galactic;
	else{
		//vlos, ma_cos(d), md => vlos, ml_cos(b), mb
		double cb = cos(b), sb = sin(b);
		double A11=(sin(decGP)*cd-cos(decGP)*sd*cos(alpha-RA_GP))/cb;
		double A12=-cos(decGP)*sin(alpha-RA_GP)/cb;
		double A21,A22;
		if(fabs(cos(lCP-l))>fabs(sin(lCP-l))){
			A21=(sd*sin(alpha-RA_GP)-sb*sin(lCP-l)*A11)/cos(lCP-l);
			A22=-(cos(alpha-RA_GP)+sb*sin(lCP-l)*A12)/cos(lCP-l);
		}else{
			A21=(cos(decGP)*cd+sin(decGP)*sd*cos(alpha-RA_GP)+sb*cos(lCP-l)*A11)/sin(lCP-l);
			A22=(sin(decGP)*sin(alpha-RA_GP)+sb*cos(lCP-l)*A12)/sin(lCP-l);
		}
	
		VecDoub GalVel {Equatorial[3],A21*Equatorial[5]+A22*Equatorial[4],
						A11*Equatorial[5]+A12*Equatorial[4]};
		for (	VecDoub::iterator it = GalVel.begin();
				it != GalVel.end(); ++it) Galactic.push_back(*it);
		return Galactic;
		}
}
std::vector<VecDoub> EquatorialToGalacticwithErrors(const VecDoub &Equatorial,
	const VecDoub &Errors){
	//alpha, dec, s => l,b,s
	double alpha = Equatorial[0], delta = Equatorial[1];
	double cd = cos(delta), sd = sin(delta);
	double b=asin(sin(decGP)*sd+cos(decGP)*cd*cos(alpha-RA_GP));
	double l=lCP-atan2(cd*sin(alpha-RA_GP),cos(decGP)*sd-sin(decGP)*cd*cos(alpha-RA_GP));
	if(l<0.)l+=2.*PI;
	if(Equatorial.size()==3){	std::vector<VecDoub> Galactic {{l,b,Equatorial[2]},Errors};
								return Galactic;}
	else{
		//vlos, ma_cos(d), md => vlos, ml_cos(b), mb
		double cb = cos(b), sb = sin(b);
		double A11=(sin(decGP)*cd-cos(decGP)*sd*cos(alpha-RA_GP))/cb;
		double A12=-cos(decGP)*sin(alpha-RA_GP)/cb;
		double A21,A22;
		if(fabs(cos(lCP-l))>fabs(sin(lCP-l))){
			A21=(sd*sin(alpha-RA_GP)-sb*sin(lCP-l)*A11)/cos(lCP-l);
			A22=-(cos(alpha-RA_GP)+sb*sin(lCP-l)*A12)/cos(lCP-l);
		}else{
			A21=(cos(decGP)*cd+sin(decGP)*sd*cos(alpha-RA_GP)+sb*cos(lCP-l)*A11)/sin(lCP-l);
			A22=(sin(decGP)*sin(alpha-RA_GP)+sb*cos(lCP-l)*A12)/sin(lCP-l);
		}
	
		std::vector<VecDoub> Galactic {
		{l,b,Equatorial[2],Equatorial[3],
		A21*Equatorial[5]+A22*Equatorial[4],A11*Equatorial[5]+A12*Equatorial[4]}
		,{Errors[0],Errors[1],Errors[2],Errors[3],
		sqrt(A21*A21*Errors[5]*Errors[5]+A22*A22*Errors[4]*Errors[4]),
		sqrt(A11*A11*Errors[5]*Errors[5]+A12*A12*Errors[4]*Errors[4])}};
		return Galactic;
		}
}
VecDoub GalacticToEquatorial(const VecDoub &Galactic){
	//l,b,s => alpha, dec, s 
	double l = Galactic[0], b = Galactic[1];
	double cb = cos(b),sb = sin(b);
	double delta=asin(cos(decGP)*cb*cos(l-lCP)+sb*sin(decGP));
	double alpha=RA_GP+atan2(cb*sin(lCP-l),sb*cos(decGP)-cb*sin(decGP)*cos(l-lCP));
	if(alpha>2.*PI)alpha-=2.*PI;
	VecDoub Equatorial {alpha,delta,Galactic[2]};
	if(Galactic.size()==3)return Equatorial;
	else{
		//vlos, ml_cos(b), mb => vlos, ma_cos(d), md
		double cd = cos(delta), sd = sin(delta);
		double A11=(sin(decGP)*cd-cos(decGP)*sd*cos(alpha-RA_GP))/cb;
		double A12=-cos(decGP)*sin(alpha-RA_GP)/cb;
		double A21,A22;
		if(fabs(cos(lCP-l))>fabs(sin(lCP-l))){
			A21=(sd*sin(alpha-RA_GP)-sb*sin(lCP-l)*A11)/cos(lCP-l);
			A22=-(cos(alpha-RA_GP)+sb*sin(lCP-l)*A12)/cos(lCP-l);
		}else{
			A21=(cos(decGP)*cd+sin(decGP)*sd*cos(alpha-RA_GP)+sb*cos(lCP-l)*A11)/sin(lCP-l);
			A22=(sin(decGP)*sin(alpha-RA_GP)+sb*cos(lCP-l)*A12)/sin(lCP-l);
		}
		double Prod = A11*A22-A12*A21;
		VecDoub EqVel {Galactic[3],(A11*Galactic[4]-A21*Galactic[5])/Prod,
						(A22*Galactic[5]-A12*Galactic[4])/Prod};
		for (	VecDoub::iterator it = EqVel.begin();
				it != EqVel.end(); ++it) Equatorial.push_back(*it);
		return Equatorial;
		}
}
}
