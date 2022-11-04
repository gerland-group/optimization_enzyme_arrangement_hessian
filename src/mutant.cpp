#include "mutant.h"

typedef Eigen::Triplet<double> T;

mutant::mutant(
	geometry* Pgeom
){
	Pg=Pgeom;
	rho.resize(Pg->Nsites);
	E=Eigen::VectorXd::Zero(Pg->Nsites);
	inv_M.resize(Pg->Nsites,Pg->Nsites);
	I.resize(Pg->Nsites, Pg->Nsites);
	I.setIdentity();
	M.resize(Pg->Nsites, Pg->Nsites);
	J=0;
}

mutant::~mutant(){
	rho.resize(0);
	E.resize(0);
	inv_M.resize(0,0);
}

mutant& mutant::operator= (
	const mutant& m
){
	rho=m.rho;
	E=m.E;
	inv_M=m.inv_M;
	J=m.J;
	return *this;
}

void mutant::rde (
){
	Eigen::SimplicialCholesky<SpM> solver;
	M = (Pg->D) - (Pg->alpha)*(SpM)E.asDiagonal();
	M.makeCompressed();
	
	solver.compute(M);
	
	inv_M = solver.solve(I);
	
	//std::ofstream oInv_M("Inv_M.out");

	//oInv_M << inv_M << std::endl;
	
	rho = inv_M*(Pg->A);
	//rho = solver.solve(Pg->A);
	for(int s=0;s<(Pg->Nsites);s++){
		if(Pg->rect_absorber.at(s)==1){
			rho(s)=0;
		}
	}
	//Calculate fluxes
	calc_react_flux();
}

void mutant::calc_react_flux (
){
	J=(Pg->dx)*(Pg->dy)*(Pg->alpha)*rho.dot(E);
}

double mutant::calc_decay_flux (
){
	double Jdecay=(Pg->dx)*(Pg->dy)*(Pg->lambda)*rho.sum();
	return Jdecay;
}

double mutant::calc_abs_flux(
){
	double Jabs=0.;
	for (int y=0; y<(Pg->Ny); ++y){
		Jabs+=(rho(y*(Pg->Nx)+(Pg->Nx)-1)+rho(y*(Pg->Nx)))/(Pg->dx)*(Pg->dy);
	}
	
	for (int x=0; x<(Pg->Nx); ++x){
		Jabs+=(rho(x)+rho(x+((Pg->Ny)-1)*(Pg->Nx)))/(Pg->dy)*(Pg->dx);
	}
	for (int s=0;s<(Pg->Nsites);s++){
		if(Pg->rect_absorber.at(s)==-1){
			Jabs+=rho(s)/(Pg->dy)*(Pg->dx);
		}else if(Pg->rect_absorber.at(s)==-2){
			Jabs+=rho(s)/(Pg->dx)*(Pg->dy);
		}
	}
	return Jabs;
}

double mutant::get_et(
){
	return E.sum()*(Pg->dx)*(Pg->dy);
}
