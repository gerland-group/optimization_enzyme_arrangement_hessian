#include <cstdlib>
#include <bits/stdc++.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "hessian.h"
#include "rectangle.h"
#include "geometry.h"
#include <eigen3/Eigen/Dense>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

//#include "sig_handler.h"

using namespace std;

bool byJ(const mutant& m1, const mutant& m2) { return (m1.J>m2.J); }

int main(){
	//Recording start and end time
	clock_t start, finish;
	start = clock();

	//Input parameters
	namespace pt = boost::property_tree;
	pt::ptree params;
	pt::read_json("params.json", params);
	
	const int Nx = params.get<int>("Nx");
	const int Ny=params.get<int>("Ny");
	const long int Nsteps=params.get<long int>("Nsteps");
	const double alpha=params.get<double>("alpha");
	const double lambda=params.get<double>("lambda");
	double det=params.get<double>("det");
	const double itol=params.get<double>("tol");
	const int Nmem=params.get<int>("Nmem");
	if(Nmem<2){
		cout << "Nmem in params.json should be at least 2" << endl;
		exit(EXIT_FAILURE);
	}	
	//Substrate sources
	vector<source> sources;
	BOOST_FOREACH (pt::ptree::value_type &src , params.get_child("sources")){
		sources.push_back(source(src));
	}
	//Absorbers
	vector<absorber> absorbers;
	boost::optional < pt::ptree& > child = params.get_child_optional( "absorbers" );
	if(child){
		BOOST_FOREACH (pt::ptree::value_type &abs , params.get_child("absorbers")){
		absorbers.push_back(absorber(abs));
		}
	}else{
		cout << "no absorber found" << endl;
	}
	
	//Derived parameters
	const double dx=1./Nx;
	const double dy=1./Ny;
	const unsigned int Nsites = Nx*Ny;
	double tol = itol;
	double time_taken = 0;
    	double ET = 0;

	//Output
	ofstream oPar("par.out");
	ofstream oEnz("E_prof.out");
	ofstream oRho("Rho_prof.out");
	//ofstream oDe("de.out");
	ofstream oTol("tol.out");
	ofstream oTrans("trans.out");
	//flux balance
	ofstream oFluxes("fluxbal.out");
	ofstream oRet("returns.out");
	ofstream oTime("time.out");
    	ofstream oComp("flux_comp.out");
	
	//Output parameters
	pt::write_json(oPar, params);	
	
	if (det<1.e-5){
		cout << "det too small for accurate hessian" << endl;
		exit(EXIT_FAILURE);
	}
	
	//loop parameters
	int x, y, //spatial coordinate
			n,    //iteration
			m;    //memory
	unsigned int s;    //lattice site
	
	//Construct geometry
	geometry geom(Nx,Ny,alpha,lambda);
	geom.build_geometry(sources,absorbers);
    	//We count the number of system sites which are not sources and not absorbers
    	int Nsites_effective = geom.count_effective_sites();
    	//We count the number of sites that have influx
    	int Nsites_with_influx = geom.count_sites_with_influx();
	
	//Solve the reaction-diffusion system with no enzymes	
	vector<mutant> before(Nmem,mutant(&geom));
	before[0].rde();
    	//Uniform is a mutant with uniform enzyme density
    	mutant uniform(&geom);
    	uniform.rde();
	//Cluster is a mutant with enzymes distributed at the sources boundary in an uniform way
	mutant cluster(&geom);
	cluster.rde();
	
	for (m=1;m<Nmem;m++){
		before[m]=before[0];
	}
	
	vector<int> pos_eq_ret;
	pos_eq_ret.reserve(Nsites);
	vector<int> sites_with_enzymes;
	sites_with_enzymes.reserve(Nsites);

	Eigen::VectorXd R = Eigen::VectorXd::Zero(Nsites);
	Eigen::VectorXd U = Eigen::VectorXd::Ones(Nsites);
	//Eigen::VectorXd rel_change = Eigen::VectorXd::Zero(Nsites);
	double * rel_change = new double [Nsites];
	double largest_with_enzymes;
	//flag for positions giving returns within tol
	bool bgradJ[Nsites];
	
	int end = 0;
	
	//We define a vector transition to keep track of the transition from a clustered enzyme configuration into an extended one.
	int transition[Nsites];
	for(s=0;s<Nsites;s++){
		transition[s]=0;
	}
	
	oTrans << "x-pos \t y-pos \t Et \t jin \t jt \t jt/jin" << endl;

	
	for(n=0; n<Nsteps; n++){
		
		//End of the loop
		if(n == Nsteps-1){
			if(end==0){
				cout << "Nsteps reached before filling the system with enzymes" << endl;
				finish = clock();
				//Time in seconds
				time_taken = double (finish-start)/double(CLOCKS_PER_SEC);
				oTime << time_taken << endl;
				exit(EXIT_FAILURE);
			}else{
				cout << "Nsteps reached and system filled with enzymes" << endl;
				finish = clock();
				//Time in seconds
				time_taken = double (finish-start)/double(CLOCKS_PER_SEC);
				oTime << time_taken << endl;
				exit(EXIT_SUCCESS);
			}
		}
		if(end == 1e4){
			cout << "System filled with enzymes" << endl;
			finish = clock();
			//Time in seconds
			time_taken = double (finish-start)/double(CLOCKS_PER_SEC);
			oTime << time_taken << endl;
			exit(EXIT_SUCCESS);
		}
		if(largest_with_enzymes+tol>1.e-1){
			cout << "tolerance too high" << endl;
			finish = clock();
			//Time in seconds
			time_taken = double (finish-start)/double(CLOCKS_PER_SEC);
	        	oTime << time_taken << endl;
			exit(EXIT_FAILURE);
		}
		//Transition data output x-pos, y-pos, Et, jin, jt, jt/jin
		//transition[s]=0 below the transition into the bulk (all enzymes are at the source at position s). If we start filling with enzymes the position adjacent to the source, we increase the value of transition[s].
		//We cannot detect transitions when sources are at the system boundaries, the adjacent position to the source would be out of bounds.
        	ET = before[0].get_et();
        	for(x=0;x<Nx;x++){
			for(y=0;y<Ny;y++){
				s=x+y*Nx;
				if((geom.rect_source.at(s)==-1)&&(transition[s]==0)){
					if(x<Nx-1){
						if((geom.rect_source.at(s+1)!=-1) &&before[0].E(s+1)>0){
							oTrans << x << " " << y << " " << before[1].get_et() << " " << -geom.A(s)*dx*dy << " " << alpha*before[1].E(s)*before[1].rho(s)*dx*dy << " " << -alpha*before[1].E(s)*before[1].rho(s)/geom.A(s) << endl;
							++transition[s];
						}
					}
					if(x>0){
						if((geom.rect_source.at(s-1)!=-1)&&before[0].E(s-1)>0){
							oTrans << x << " " << y << " " << before[1].get_et() << " " << -geom.A(s)*dx*dy << " " << alpha*before[1].E(s)*before[1].rho(s)*dx*dy << " " << -alpha*before[1].E(s)*before[1].rho(s)/geom.A(s) << endl;
							++transition[s];
						}
					}
					if(y>0){
						if((geom.rect_source.at(s-Nx)!=-1)&&before[0].E(s-Nx)>0){
							oTrans << x << " " << y << " " << before[1].get_et() << " " << -geom.A(s)*dx*dy << " " << alpha*before[1].E(s)*before[1].rho(s)*dx*dy << " " << -alpha*before[1].E(s)*before[1].rho(s)/geom.A(s) << endl;
							++transition[s];
						}
					}
					if(y<Ny-1){
						if((geom.rect_source.at(s+Nx)!=-1)&&before[0].E(s+Nx)>0){
							oTrans << x << " " << y << " " << before[1].get_et() << " " << -geom.A(s)*dx*dy << " " << alpha*before[1].E(s)*before[1].rho(s)*dx*dy << " " << -alpha*before[1].E(s)*before[1].rho(s)/geom.A(s) << endl;
							++transition[s];
						}
					}
				}
			}
		}
		
		//Output the enzyme, the rho profiles and the returns in matrix form for approx 1000 times. Top line in the output is the bottom line in the system, the top line in the system is the bottom line in the output (gnuplot matrix format)
		//We also monitor the reaction flux at the boundary, in the positions where there is influx.
		//if(n%int(Nsteps/1000)==0){
		if(n%6000==0){
			for(y=0;y<Ny;y++){				
				for(x=0;x<Nx;x++){
					s=x+y*Nx;
					oEnz << before[0].E(s) << " ";
					oRho << before[0].rho(s) << " ";
					oRet << R(s) << " ";
                    			//The uniform mutant has a uniform enzyme density among the Nsites_effective
					if(geom.rect_source.at(s)!=1 || geom.rect_absorber.at(s)!=1){
						uniform.E(s) = ET*Nx*Ny/Nsites_effective;
					}
					//The cluster mutant has enzymes only near the source boundaries
					if(geom.rect_source.at(s)==-1){
						cluster.E(s) = ET*Nx*Ny/Nsites_with_influx;
					}
				}
				oEnz << endl;
				oRho << endl;
				oRet << endl;
			}
			oEnz << endl;
			oRho << endl;
			oRet << endl;
			
			//Solve for a uniform enzyme density
			uniform.rde();
			//Solve for the cluster mutant
			cluster.rde();
			//Compare uniform, cluster and optimal fluxes
			oComp << ET << " " << before[0].J << " " << uniform.J << " " << cluster.J << endl;
			//Check flux balance
			oFluxes << ET << " " << geom.calc_influx() << " " << before[0].J << " " << before[0].calc_abs_flux() << " " << before[0].calc_abs_flux()+before[0].J << endl;
			//tol accepted as we increase et
			oTol << ET << " " << largest_with_enzymes+tol << endl;
		}
		
		//If system is full we should update all sites, so skip site selection 
		if(sites_with_enzymes.size()<Nsites) {
			//Resetting tol to itol
			tol = itol;
			
			//Determine the returns vector and its maximum
			R = (U + alpha*before[0].inv_M.transpose()*before[0].E).array()*before[0].rho.array()*alpha*dx*dy;
			Eigen::VectorXd::Index max_ind;
			double Rmax = R.maxCoeff(&max_ind);
			
			for(s=0;s<Nsites;s++){
				rel_change[s] = fabs(1.-R(s)/Rmax);
			}
			//rel_change = (U - R/Rmax).cwiseAbs() ;
	
			largest_with_enzymes=0.;
			pos_eq_ret.clear();
			pos_eq_ret=sites_with_enzymes;
			for(s=0; s<Nsites; ++s) {
				bgradJ[s]=0;
			}
			//Find out is the largest deviation in returns among sites that contain enzymes
			if(sites_with_enzymes.size()>0) {
				for(unsigned int i=0; i<sites_with_enzymes.size(); ++i){
					s=sites_with_enzymes.at(i);
					if(rel_change[s]>largest_with_enzymes) {
						largest_with_enzymes=rel_change[s];
					}
					bgradJ[s] = 1 ;
				}
			}
			//After we have included all sites with enzymes, add any other sites within tolerance
			for(s=0; s<Nsites; ++s){
				if(!bgradJ[s]) {
					if(rel_change[s] < largest_with_enzymes+tol){
						pos_eq_ret.push_back(s);
						sites_with_enzymes.push_back(s);
					}
				}
			}
			sort(sites_with_enzymes.begin(),sites_with_enzymes.end());
			sort(pos_eq_ret.begin(),pos_eq_ret.end());
		}
		
		//Set the new starting system and set the mutants before
		m=Nmem-1;
		while(m!=0){
			before[m] = before[m-1];
			--m;
		}
		if (pos_eq_ret.size()==1){
			//All det enzymes at one position
			before[0].E(pos_eq_ret.at(0)) += det;
		}else{
			if(pos_eq_ret.size() == Nsites){
				++end;
			}
		
			//Compute the hessian 
			Eigen::MatrixXd H = Eigen::MatrixXd(pos_eq_ret.size(),pos_eq_ret.size());
			hessian(H, U, before[0], alpha, Nsites, pos_eq_ret);
			
			//Solve the system (1) = H (de/det)
			//Euler update
			Eigen::VectorXd de(pos_eq_ret.size());
			
			solve_de(de,H);
			
/*			//Runge-Kutta 4
			//Set the systems needed for RK4
			Eigen::VectorXd * RK = new Eigen::VectorXd[4];
			mutant step;
			int index_RK;
			for(index_RK=0;index_RK<4;index_RK++){
				RK[index_RK].resize(count_eq_returns); 
			}
			
			solve_de(RK[0],H);
			
			for(index_RK=1;index_RK<4;index_RK++){
				step = before[0];
				for(index_eq_ret=0;index_eq_ret<count_eq_returns;index_eq_ret++){
					step.E(pos_eq_ret[index_eq_ret])+= RK[index_RK-1](index_eq_ret)*det/((index_RK==3) ? 1.:2.);
				}
				step.rde();
				
				hessian(H, U, count_eq_returns, step, alpha, Nx, Ny, pos_eq_ret);
				solve_de(RK[index_RK],H);
			}
*/
			
			//Euler update
			for(unsigned int i=0; i<pos_eq_ret.size(); ++i){
				//if(n%1000==0){
					//oDe << pos_eq_ret.at(i) << " " << de(i) << " ";
				//}
				before[0].E(pos_eq_ret.at(i))+=de(i)*det;
			}
			//oDe << endl;
			
/*			//Runge-Kutta 4
			for(index_RK=0;index_RK<count_eq_returns;index_RK++){
				oDe << 1./6.*det*(RK[0](index_RK)+2.*RK[1](index_RK)+2.*RK[2](index_RK)+RK[3](index_RK)) << " ";
				before[0].E(pos_eq_ret[index_RK])+= 1./6.*det*(RK[0](index_RK)+2.*RK[1](index_RK)+2.*RK[2](index_RK)+RK[3](index_RK));
			} oDe << endl;
			before[0].rde();
*/			
			
			H.resize(0,0);
			de.resize(0);
			//for(index_RK=0;index_RK<4;index_RK++){
				//RK[index_RK].resize(0);
			//}
			//delete[] RK;
			// delete[] pos_eq_ret;
		}
		before[0].rde();
	}
	
	//Free memory
	R.resize(0);
	U.resize(0);
	delete[] rel_change;
	//rel_change.resize(0);
	exit(EXIT_SUCCESS);
}
