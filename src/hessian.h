#ifndef HESSIAN_H
#define HESSIAN_H

#include <iostream>
#include <vector>
#include "mutant.h"
#include <eigen3/Eigen/Dense>

void hessian (Eigen::MatrixXd&, Eigen::VectorXd&, const mutant& , const double, const int, const std::vector<int>& );
void solve_de (Eigen::VectorXd&, Eigen::MatrixXd&);

#endif
