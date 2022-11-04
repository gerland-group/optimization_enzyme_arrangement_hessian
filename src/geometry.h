#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <eigen3/Eigen/Sparse>
#include <vector>
#include "rectangle.h"

typedef Eigen::SparseMatrix<double> SpM;

class geometry{
	public:
	Eigen::VectorXd A;
	SpM D;
	int Nx;
	int Ny;
	int Nsites;
	double dx;
	double dy;
//model parameters
	double alpha;
	double lambda;
//rect_source[i]=1 if the i-site is within the source; rect_source[i]=-1 if there is influx at the i-site and =0 otherwise
//rect_absorber[i]=1 if the site i is part of the absorber, rect_absorber[i]=-1 in the surrounding sites below and above the absorber, =-2 on the left and right of the absorber, and =0 otherwise	
	std::vector<int> rect_source;
	std::vector<int> rect_absorber;
	
	geometry(int, int, double, double);
	~geometry();

	void build_geometry(const std::vector<source>&, const std::vector<absorber>&);
	template <class Type> void check_bounds(const std::vector<Type>&);
	void fill_D_matrix();
	double calc_influx();

	void set_influx(const source&);
	void set_absorber(const absorber&);
    int count_effective_sites();
    int count_sites_with_influx();
};

#endif
