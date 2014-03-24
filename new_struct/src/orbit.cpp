/*======================================*/
/* 				  Orbits	 			*/
/*======================================*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface.h"
#include "utils.h"
#include "gnuplot/gnuplot_i.h"
#include "potential.h"
#include "orbit.h"

int derivs(double t, const double y[], double f[],void *params){
	/* derivatives for orbit integration */
	Potential_JS *P = (Potential_JS*) params;
	VecDoub x = {y[0],y[1],y[2]};
	VecDoub F = P->Forces(x);
	for(int i=0;i<3;i++){ f[i]=y[i+3]; f[i+3]=F[i]; }
	return GSL_SUCCESS;
}

Orbit::Orbit(Potential_JS *P, double eps){
	/* Create instance of orbit with associated Potential_JS P */
	O = new ode(derivs,6,eps, P);
	labels = {"x","y","z","vx","vy","vz"};
}

VecDoub Orbit::integrate(VecDoub x, double t_interval, double step_size){
	/* Integrate from time t=0 to t=t_interval at regular steps step_size */
	/* and store results in the array results							  */
	double y[6]; VecDoub xf(6); for(int i=0;i<6;i++)y[i]=x[i];
	double start=0.;
	results.clear();
	results.push_back(x);
	while((t_interval-start)>0.){
		O->step(start,start+step_size,y,step_size);
		for(int i=0;i<6;i++) xf[i]=y[i];
		results.push_back(xf);
		start+=step_size;
	}
	return xf;
}
// void Orbit::plot(int i, int j, std::string name){
// 	/* Plot the results of the orbit integration  */
// 	/* i,j give the indices of the arrays to plot */
// 	VecDoub x,y; std::string x_s=labels[i], y_s=labels[j];
// 	for(unsigned int k=0;k<results.size();k++){
// 		x.push_back(results[k][i]);y.push_back(results[k][j]);
// 	} 
// 	Gnuplot G("lines ls 1");
// 	G.set_xrange(1.1*Min(x),1.1*Max(x)).set_yrange(1.1*Min(y),1.1*Max(y));
// 	G.set_xlabel(x_s).set_ylabel(y_s);
// 	G.savetotex(name).plot_xy(x,y);
// 	G.outputpdf(name);
// }
