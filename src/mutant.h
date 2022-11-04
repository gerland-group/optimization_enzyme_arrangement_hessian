#ifndef MUTANT_H
#define MUTANT_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <cmath>
#include "geometry.h"

typedef Eigen::SparseMatrix<double> SpM;

class mutant{
	private:
		geometry* Pg;
		SpM M;
		Eigen::MatrixXd I;
	public:
		Eigen::VectorXd E;
		Eigen::VectorXd rho;
		Eigen::MatrixXd inv_M;
	double J;

	mutant(geometry*);
	~mutant();
	mutant& operator=(const mutant&);

	void rde ();
	void calc_react_flux();
	double calc_decay_flux ();
	double calc_abs_flux();
	double get_et();
};

#endif
