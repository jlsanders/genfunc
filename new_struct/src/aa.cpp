/*======================================*/
/*Actions, Angles, Frequencies & Hessian*/
/*======================================*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "GSLInterface.h"
#include "utils.h"
#include "potential.h"
#include "aa.h"

// ==========================================================================================
// Prolate Stackel Angle-action calculator
// ==========================================================================================

double ptau2ROOT_AxiSym(double tau, void *params){
	/* for finding roots of p_tau^2*2.0*(tau+Alpha)  */
	root_struct_axi *RS = (root_struct_axi *) params;
	return (RS->Ints[0]-RS->Ints[1]/(tau+RS->P->alpha())-RS->Ints[2]/(tau+RS->P->gamma())+RS->P->G(tau));
	}

VecDoub  Actions_AxisymmetricStackel::find_limits(VecDoub tau, VecDoub ints){

	double lambda = tau[0], nu = tau[2];
	VecDoub limits;
	// create a structure to store parameters for ptau2ROOT
	root_struct_axi RS(Pot,ints);
	// find roots of p^2(lam)
	double laminner=lambda;
	while(ptau2ROOT_AxiSym(laminner, &RS)>0.0)	laminner-=.1*(laminner+Pot->alpha());
	double lamouter=lambda;
	while(ptau2ROOT_AxiSym(lamouter, &RS)>0.)	lamouter*=1.1;

	root_find RF(1e-5,100);
	limits.push_back(RF.findroot(&ptau2ROOT_AxiSym,laminner,lambda,&RS));
	limits.push_back(RF.findroot(&ptau2ROOT_AxiSym,lambda,lamouter,&RS));
	limits.push_back(-Pot->gamma()+TINY);
	// find root of p^2(nu)
	double nuouter=nu;
	while(ptau2ROOT_AxiSym(nuouter, &RS)<0.)	nuouter+=0.1*(-Pot->alpha()-nuouter);
	limits.push_back(RF.findroot(&ptau2ROOT_AxiSym,nu,nuouter,&RS));
	return limits;
}

double ptau2_AxiSym(double tau, void *params){
	//p^2(tau) using global integrals
	root_struct_axi *RS = (root_struct_axi *) params;
	return (RS->Ints[0]-RS->Ints[1]/(tau+RS->P->alpha())-RS->Ints[2]/(tau+RS->P->gamma())+RS->P->G(tau))
			/(2.0*(tau+RS->P->alpha()));
}
	
double J_integrand_AxiSym(double theta, void *params){
	// need to set taubargl and Deltagl
	action_struct_axi * AS = (action_struct_axi *) params;
	double tau=AS->taubargl+AS->Deltagl*sin(theta);
	return sqrt(MAX(0.,ptau2_AxiSym(tau,&(AS->RS))))*cos(theta);
}

VecDoub Actions_AxisymmetricStackel::actions(VecDoub x){
	VecDoub tau = Pot->xv2tau(x);
	VecDoub integrals = Pot->x2ints(x,&tau);
	VecDoub limits = find_limits(tau,integrals);
	VecDoub actions;
	// JR
	double taubar = 0.5*(limits[0]+limits[1]);
	double Delta = 0.5*(limits[1]-limits[0]);
	action_struct_axi AS(Pot,integrals,taubar,Delta);
	actions.push_back(Delta*GaussLegendreQuad(&J_integrand_AxiSym,-.5*PI,.5*PI,&AS)/PI);
	// Lz
	actions.push_back(sqrt(2.0*integrals[1]));
	// Jz
	taubar = 0.5*(limits[2]+limits[3]);
	Delta = 0.5*(limits[3]-limits[2]);
	AS = action_struct_axi (Pot,integrals,taubar,Delta);
	actions.push_back(2.*Delta*GaussLegendreQuad(&J_integrand_AxiSym,-.5*PI,.5*PI,&AS)/PI);
	return actions;
}

// ==========================================================================================
// Triaxial Stackel Angle-action calculator
// ==========================================================================================

double ptau2ROOT_Triax(double tau, void *params){
	/* for finding roots of p_tau^2*2.0*(tau+Alpha)  */
	root_struct_triax *RS = (root_struct_triax *) params;
	return (RS->Ints[0]-RS->Ints[1]/(tau+RS->P->alpha())-RS->Ints[2]/(tau+RS->P->gamma())+RS->P->G(tau))
			/(2.0*(tau+RS->P->beta()));
	}

VecDoub Actions_TriaxialStackel::find_limits(VecDoub tau,VecDoub ints){

	double lambda = tau[0], mu=tau[1], nu = tau[2];
	root_find RF(1e-5,100);
	VecDoub limits;
	// create a structure to store parameters for ptau2ROOT
	root_struct_triax RS(Pot,ints);

	// find roots of p^2(lam)
	double laminner=lambda, lamouter=lambda;
	while(ptau2ROOT_Triax(laminner, &RS)>0.0 and (laminner+Pot->alpha())>1e-5)	laminner-=.1*(laminner+Pot->alpha());
	if((laminner+Pot->alpha())>1e-5) limits.push_back(RF.findroot(&ptau2ROOT_Triax,laminner,lambda,&RS));
	else limits.push_back(-Pot->alpha());
	while(ptau2ROOT_Triax(lamouter, &RS)>0.)	lamouter*=1.1;
	limits.push_back(RF.findroot(&ptau2ROOT_Triax,lambda,lamouter,&RS));

	// find root of p^2(mu)
	double muinner=mu, muouter=mu;
	while(ptau2ROOT_Triax(muinner, &RS)>0. and (muinner+Pot->beta())>1e-5)	muinner-=.1*(muinner+Pot->beta());
	if((muinner+Pot->beta())>1e-5) limits.push_back(RF.findroot(&ptau2ROOT_Triax,muinner,mu,&RS));
	else limits.push_back(-Pot->beta());
	while(ptau2ROOT_Triax(muouter, &RS)>0. and (muouter+Pot->alpha())<-1e-5)	muouter+=0.1*(-Pot->alpha()-muouter);
	if((muouter+Pot->alpha())<-1e-5) limits.push_back(RF.findroot(&ptau2ROOT_Triax,mu,muouter,&RS));
	else limits.push_back(-Pot->alpha());

	// find root of p^2(nu)
	double nuinner=nu, nuouter=nu;
	while(ptau2ROOT_Triax(nuinner, &RS)>0. and (nuinner+Pot->gamma())>1e-5)	nuinner-=.1*(nuinner+Pot->gamma()); 
	if((nuinner+Pot->gamma())>1e-5) limits.push_back(RF.findroot(&ptau2ROOT_Triax,nuinner,nu,&RS));
	else limits.push_back(-Pot->gamma());
	while(ptau2ROOT_Triax(nuouter, &RS)>0. and (nuouter+Pot->beta())<-1e-5)	nuouter+=0.1*(-Pot->beta()-nuouter);
	if((nuouter+Pot->beta())<-1e-5) limits.push_back(RF.findroot(&ptau2ROOT_Triax,nu,nuouter,&RS));
	else limits.push_back(-Pot->beta());

	return limits;
}
	
double J_integrand_Triax(double theta, void *params){
	// need to set taubargl and Deltagl
	action_struct_triax * AS = (action_struct_triax *) params;
	double tau=AS->taubargl+AS->Deltagl*sin(theta);
	double ptau = ((tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->Ints[0]
			-AS->Ints[1]*(tau+AS->P->gamma())-AS->Ints[2]*(tau+AS->P->alpha())
			+(tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->P->G(tau))
			/(2.0*(tau+AS->P->alpha())*(tau+AS->P->beta())*(tau+AS->P->gamma()));
	return sqrt(MAX(0.,ptau))*cos(theta);
}

double dJdH_integrand_Triax(double theta, void *params){
	// need to set taubargl and Deltagl
	action_struct_triax * AS = (action_struct_triax *) params;
	double tau=AS->taubargl+AS->Deltagl*sin(theta);
	double ptau = ((tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->Ints[0]
			-AS->Ints[1]*(tau+AS->P->gamma())-AS->Ints[2]*(tau+AS->P->alpha())
			+(tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->P->G(tau))
			/(2.0*(tau+AS->P->alpha())*(tau+AS->P->beta())*(tau+AS->P->gamma()));
	return sqrt(MAX(0.,1./ptau))*cos(theta)/(tau+AS->P->beta());
}

double dJdI2_integrand_Triax(double theta, void *params){
	// need to set taubargl and Deltagl
	action_struct_triax * AS = (action_struct_triax *) params;
	double tau=AS->taubargl+AS->Deltagl*sin(theta);
	double ptau = ((tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->Ints[0]
			-AS->Ints[1]*(tau+AS->P->gamma())-AS->Ints[2]*(tau+AS->P->alpha())
			+(tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->P->G(tau))
			/(2.0*(tau+AS->P->alpha())*(tau+AS->P->beta())*(tau+AS->P->gamma()));
	return -sqrt(MAX(0.,1./ptau))*cos(theta)/(tau+AS->P->beta())/(tau+AS->P->alpha());
}

double dJdI3_integrand_Triax(double theta, void *params){
	// need to set taubargl and Deltagl
	action_struct_triax * AS = (action_struct_triax *) params;
	double tau=AS->taubargl+AS->Deltagl*sin(theta);
	double ptau = ((tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->Ints[0]
			-AS->Ints[1]*(tau+AS->P->gamma())-AS->Ints[2]*(tau+AS->P->alpha())
			+(tau+AS->P->alpha())*(tau+AS->P->gamma())*AS->P->G(tau))
			/(2.0*(tau+AS->P->alpha())*(tau+AS->P->beta())*(tau+AS->P->gamma()));
	return -sqrt(MAX(0.,1./ptau))*cos(theta)/(tau+AS->P->beta())/(tau+AS->P->gamma());
}

VecDoub Actions_TriaxialStackel::actions(VecDoub x, bool freq_yes){
	VecDoub tau = Pot->xv2tau(x);
	VecDoub integrals = Pot->tau2ints(tau);
	VecDoub limits = find_limits(tau,integrals);
	VecDoub actions, freqs;

	// We need to check which coordinates are oscillating and which are circulating
	// and multiply by the appropriate factor
	VecDoub circ={1.,1.,1.};
	if(limits[0]==-Pot->alpha()) circ[0] = 1.; else circ[0] = 0.5;
	// if(limits[2]==-Pot->beta() and limits[3]==-Pot->alpha()) circ[1] = 1.; else circ[1] = 0.5;
	// if(limits[4]==-Pot->gamma() and limits[5]==-Pot->beta()) circ[2] = 0.5; else circ[2] = 1.;

	// JR
	double taubar = 0.5*(limits[0]+limits[1]);
	double Delta = 0.5*(limits[1]-limits[0]);
	action_struct_triax AS(Pot,integrals,taubar,Delta);
	actions.push_back(2.*circ[0]*Delta*GaussLegendreQuad(&J_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);

	if(freq_yes){
		freqs.push_back(0.5*circ[0]*Delta*GaussLegendreQuad(&dJdH_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
		freqs.push_back(0.5*circ[0]*Delta*GaussLegendreQuad(&dJdI2_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
		freqs.push_back(0.5*circ[0]*Delta*GaussLegendreQuad(&dJdI3_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
	}

	// Jp
	taubar = 0.5*(limits[2]+limits[3]);
	Delta = 0.5*(limits[3]-limits[2]);
	AS = action_struct_triax(Pot,integrals,taubar,Delta);
	actions.push_back(2.*circ[1]*Delta*GaussLegendreQuad(&J_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);

	if(freq_yes){
		freqs.push_back(0.5*circ[1]*Delta*GaussLegendreQuad(&dJdH_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
		freqs.push_back(0.5*circ[1]*Delta*GaussLegendreQuad(&dJdI2_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
		freqs.push_back(0.5*circ[1]*Delta*GaussLegendreQuad(&dJdI3_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
	}

	// Jz
	taubar = 0.5*(limits[4]+limits[5]);
	Delta = 0.5*(limits[5]-limits[4]);
	AS = action_struct_triax(Pot,integrals,taubar,Delta);
	actions.push_back(2.*circ[2]*Delta*GaussLegendreQuad(&J_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);

	if(freq_yes){
		freqs.push_back(0.5*circ[2]*Delta*GaussLegendreQuad(&dJdH_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
		freqs.push_back(0.5*circ[2]*Delta*GaussLegendreQuad(&dJdI2_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
		freqs.push_back(0.5*circ[2]*Delta*GaussLegendreQuad(&dJdI3_integrand_Triax,-.5*PI,.5*PI,&AS)/PI);
		double det = freqs[0]*(freqs[4]*freqs[8]-freqs[5]*freqs[7])
					-freqs[1]*(freqs[3]*freqs[8]-freqs[5]*freqs[6])
					+freqs[2]*(freqs[3]*freqs[7]-freqs[4]*freqs[6]);
		actions.push_back((freqs[4]*freqs[8]-freqs[5]*freqs[7])/det);
		actions.push_back((freqs[7]*freqs[2]-freqs[8]*freqs[1])/det);
		actions.push_back((freqs[1]*freqs[5]-freqs[2]*freqs[4])/det);
	}

	return actions;
}
// ==========================================================================================
// aa.cpp
