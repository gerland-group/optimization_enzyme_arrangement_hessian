#ifndef RECTANGLE_H
#define RECTANGLE_H

#include <vector>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
class geometry;

class rectangle{
	public:
	int x_pos;
	int y_pos;
	int x_size;
	int y_size;
	rectangle(const boost::property_tree::ptree::value_type&);
};

class source: public rectangle{
	public:
	double A_0;
	
	//~source();
	source(const boost::property_tree::ptree::value_type&);
};

class absorber: public rectangle{
	public:
	absorber(const boost::property_tree::ptree::value_type&);
};

#endif
