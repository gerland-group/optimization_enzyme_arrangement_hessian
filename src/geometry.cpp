#include "geometry.h"

typedef Eigen::Triplet<double> T;

geometry::geometry(
	int N_sites_x,
	int N_sites_y,
	double reaction_rate,
	double decay_rate
){
	Nx=N_sites_x;
	Ny=N_sites_y;
	alpha=reaction_rate;
	lambda=decay_rate;
	//Derived parameters
	dx=1./Nx;
	dy=1./Ny;
	Nsites=Nx*Ny;
	
	D.resize(Nsites,Nsites);
	D.reserve(Eigen::VectorXd::Constant(Nsites,5));
	A=Eigen::VectorXd::Zero(Nsites);
	rect_source.resize(Nsites,0);
	rect_absorber.resize(Nsites,0);
}

geometry::~geometry(){
	A.resize(0);
	D.resize(0,0);
	rect_source.resize(0);
	rect_absorber.resize(0);
}

void geometry::build_geometry(
	const std::vector<source>& sources,
	const std::vector<absorber>& absorbers
){
 	
	int x,y;
	
	std::ofstream oA ("A.out");
	//std::ofstream oD ("D.out");
	//Influx vector A: rectangular sources with linear influx with intensity A_0 on each side of the rectangles
	//Check sources bounds
	check_bounds(sources);
	for(std::vector<source>::size_type i = 0; i != sources.size(); i++) {
		set_influx(sources[i]);
	}
	A/=(dx*dy);
	
	//Check absorbers bounds
	check_bounds(absorbers);
	//Setting the absorbers
	for(std::vector<absorber>::size_type i = 0; i != absorbers.size(); i++) {
		set_absorber(absorbers[i]);
	}
	
	//Check for overlaps between absorbers and sources
	for(std::vector<int>::size_type s=0;s != rect_source.size(); s++){
		if(rect_source.at(s)==1 && rect_absorber.at(s)==1){
			std::cout << "there is an overlap between absorbers and sources at position (x,y)=" << (s%Nx) << "," << (s/Nx) << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
	
	//Build The D matrix
	fill_D_matrix();

	//Output D
	//std::cout << D << std::endl;

	//Output vector A
	for(y=0;y<Ny;y++){
		for(x=0;x<Nx;x++){
			oA << A(x+Nx*y) << " " ;
		}
		oA << std::endl;
	}
	oA << std::endl;

	std::cout << "geometrical properties set" << std::endl;
}

template <class Type>
void geometry::check_bounds(
	const std::vector<Type>& Types
){
	for(typename std::vector<Type>::size_type i = 0; i != Types.size(); i++){
		if(Types[i].x_size > Nx || Types[i].y_size > Ny || Types[i].x_size < 0 || Types[i].y_size < 0){
			std::cout << "choose a valid size for the rectangular " << typeid(Types[i]).name() << " " << std::endl;
			std::cout << "x_size= " << Types[i].x_size << " and y_size= " << Types[i].y_size << " should be in the intervals [" << 0 << "," << Nx << "] and [" << 0 << "," << Ny << "] respectively" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		if((Types[i].x_pos > (Nx-Types[i].x_size))||(Types[i].x_pos < 0)){
			std::cout << "choose a valid x-position for the rectangular " << typeid(Types[i]).name() << " " << std::endl;
			std::cout << "x_pos= " << Types[i].x_pos << " should be in the interval [" << 0 << "," << Nx-Types[i].x_size << "]" << std::endl;
			std::exit(EXIT_FAILURE);
		}
		
		if((Types[i].y_pos > (Ny-Types[i].y_size))||(Types[i].y_pos < 0)){
			std::cout << "choose a valid y-position for the rectangular " << typeid(Types[i]).name() << " " << std::endl;
			std::cout << "y_pos= " << Types[i].y_pos << " should be in the interval [" << 0 << "," << Ny-Types[i].y_size << "]" << std::endl;
			std::exit(EXIT_FAILURE);
		}
	}
}

void geometry::fill_D_matrix(){
	
	std::vector<T> tripletList;
	tripletList.reserve(Nsites+2*(Nsites-1)+2*(Ny));
	
	int x,y,s;
	
	for(y=0;y<Ny;y++){
		for(x=0;x<Nx;x++){
			s=x+y*Nx;
			if(rect_source.at(s)==1 || rect_absorber.at(s)==1){
				tripletList.push_back(T(s,s,1.));
				A[s]=0;
			}else{
				if(x<Nx-1 && rect_source.at(s+1)!=1 && rect_absorber.at(s+1)!=1){
					tripletList.push_back(T(s,s+1,1./pow(dx,2)));
				}
				if(x>0 && rect_source.at(s-1)!=1 && rect_absorber.at(s-1)!=1){
					tripletList.push_back(T(s,s-1,1./pow(dx,2)));
				}
				if(y<Ny-1 && rect_source.at(s+Nx)!=1 && rect_absorber.at(s+Nx)!=1){
					tripletList.push_back(T(s,s+Nx,1./pow(dy,2)));
				}
				if(y>0 && rect_source.at(s-Nx)!=1 && rect_absorber.at(s-Nx)!=1){
					tripletList.push_back(T(s,s-Nx,1./pow(dy,2)));
				}
				//With reflecting boundary conditions on the outer boundary
				//tripletList.push_back(T(s,s,-(2-(x==0?1:0)-(x==Nx-1?1:0))/pow(dx,2)-(2-(y==0?1:0)-(y==Ny-1?1:0))/pow(dy,2)-lambda));
				//std::cout << "reflective boundary conditions on the outer boundary" << std::endl;
				//With absorbing boundary conditions on the outer boundary
				tripletList.push_back(T(s,s,-(2-(((x<Nx-1)&&rect_source.at(s+1)==1)?1:0)-(((x>0)&&rect_source.at(s-1)==1)?1:0))/pow(dx,2)- (2-(((y<Ny-1)&&rect_source.at(s+Nx)==1)?1:0)-(((y>0)&&rect_source.at(s-Nx)==1)?1:0))/pow(dy,2)-lambda));
				//std::cout << "absorbing boundary conditions on the outer boundary" << std::endl;
				//With absorbing left and right boundaries and reflecting at the top and bottom
				//tripletList.push_back(T(s,s,-(2-(((x<Nx-1)&&rect_source.at(s+1)==1)?1:0)-(((x>0)&&rect_source.at(s-1)==1)?1:0))/pow(dx,2)-(2-((y==0)?1:0)-((y==Ny-1)?1:0)-(((y<Ny-1)&&rect_source.at(s+Nx)==1)?1:0)-(((y>0)&&rect_source.at(s-Nx)==1)?1:0))/pow(dy,2)-lambda));
				//std::cout << "absorbing: right, left; reflecting: tob, bottom" << std::endl;
				//(Beta) Influx along external left boundary, reflecting top and bottom, absorbing right
				//tripletList.push_back(T(s,s,-(((x==0)?1:2)-(((x<Nx-1)&&rect_source.at(s+1)==1)?1:0)-(((x>0)&&rect_source.at(s-1)==1)?1:0))/pow(dx,2)-(2-((y==0)?1:0)-((y==Ny-1)?1:0)-(((y<Ny-1)&&rect_source.at(s+Nx)==1)?1:0)-(((y>0)&&rect_source.at(s-Nx)==1)?1:0))/pow(dy,2)-lambda));
				//std::cout << "influx: left; reflecting: top, bottom; absorbing: right" << std::endl;
			}
		}
	}
	
	D.setFromTriplets(tripletList.begin(),tripletList.end());
	D.makeCompressed();
	
	//oD << D << std::endl;
}

double geometry::calc_influx(
){
	double Jin=0.;
	for(int s=0;s<Nsites;s++){
		Jin+=(-A(s)*dx*dy);
	}
	
	return Jin;
}

void geometry::set_influx (
	const source& s
){
	
	int x,y; //Spatial coordinates relative to the x/y position.
	
	int or_pos_x, or_pos_y,
		or_pos_lattice;
	
	for(x=0;x<s.x_size;x++){
		for(y=0;y<s.y_size;y++){
			or_pos_x = s.x_pos + x;
			or_pos_y = s.y_pos + y;
			or_pos_lattice = or_pos_x + Nx*or_pos_y;
			
			rect_source.at(or_pos_lattice)=1;
			
			if(y==0 && (or_pos_y-1)>=0 && rect_source.at(or_pos_x+ Nx*(or_pos_y-1))!=1){
				A(or_pos_x+ Nx*(or_pos_y-1))-=(s.A_0);
				rect_source.at(or_pos_x+ Nx*(or_pos_y-1))=-1;
			}
			if(y==s.y_size-1 && (or_pos_y+1)<Ny && rect_source.at(or_pos_x+ Nx*(or_pos_y+1))!=1){
				A(or_pos_x+ Nx*(or_pos_y+1))-=(s.A_0);
				rect_source.at(or_pos_x+ Nx*(or_pos_y+1))=-1;
			}
			if(x==0 && (or_pos_x-1)>=0 && rect_source.at(or_pos_x-1 + Nx*or_pos_y)!=1){
				A(or_pos_x-1 + Nx*or_pos_y)-=(s.A_0);
				rect_source.at(or_pos_x-1 + Nx*or_pos_y)=-1;
			}
			if(x==s.x_size-1 && (or_pos_x+1)<Nx && rect_source.at(or_pos_x+1 + Nx*or_pos_y)!=1){
				A(or_pos_x+1 + Nx*or_pos_y)-=(s.A_0);
				rect_source.at(or_pos_x+1 + Nx*or_pos_y)=-1;
			}
		}
	}
}


void geometry::set_absorber(
	const absorber& a
) {
	
	int x,y; //Spatial coordinates relative to the x/y position.
	
	int or_pos_x, or_pos_y,
		or_pos_lattice;
		
	for(x=0;x<a.x_size;x++){
		for(y=0;y<a.y_size;y++){
			or_pos_x = a.x_pos + x;
			or_pos_y = a.y_pos + y;
			or_pos_lattice = or_pos_x + Nx*or_pos_y;
			
			rect_absorber.at(or_pos_lattice)=1;
			if(y==0 && (or_pos_y-1)>=0 && rect_absorber.at(or_pos_x+ Nx*(or_pos_y-1))!=1){
				rect_absorber.at(or_pos_x+ Nx*(or_pos_y-1))=-1;
			}
			if(y==a.y_size-1 && (or_pos_y+1)<Ny && rect_absorber.at(or_pos_x+ Nx*(or_pos_y+1))!=1){
				rect_absorber.at(or_pos_x+ Nx*(or_pos_y+1))=-1;
			}
			if(x==0 && (or_pos_x-1)>=0 && rect_absorber.at(or_pos_x-1 + Nx*or_pos_y)!=1){
				rect_absorber.at(or_pos_x-1 + Nx*or_pos_y)=-2;
			}
			if(x==a.x_size-1 && (or_pos_x+1)<Nx && rect_absorber.at(or_pos_x+1 + Nx*or_pos_y)!=1){
				rect_absorber.at(or_pos_x+1 + Nx*or_pos_y)=-2;
			}
		}
	}
}

int geometry::count_effective_sites(
){
    int Nsites_effective=0;

	for(int s=0;s<Nsites;s++){
		if(rect_source.at(s)!=1 || rect_absorber.at(s)!=1){
			Nsites_effective++;
        }
    }

    return Nsites_effective;
}

int geometry::count_sites_with_influx(
){
	int Nsites_with_influx=0;
	for(int s=0;s<Nsites;s++){
		if(rect_source.at(s)==-1){
			Nsites_with_influx++;
		}
	}
	
	return Nsites_with_influx;
}

