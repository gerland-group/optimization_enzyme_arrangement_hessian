#include "rectangle.h"
#include "geometry.h"

rectangle::rectangle(
	const boost::property_tree::ptree::value_type& tree
){
	x_pos=tree.second.get<int>("x_pos");
	y_pos=tree.second.get<int>("y_pos");
	x_size=tree.second.get<int>("x_size");
	y_size=tree.second.get<int>("y_size");
}

source::source(
	const boost::property_tree::ptree::value_type& sources_tree
):rectangle(sources_tree){
	A_0=sources_tree.second.get<double>("A_0");
	
	if(A_0 <= 0){
		std::cout << "choose a valid influx strength" << std::endl;
		std::cout << "A_0= " << A_0 << " should be positive" << std::endl;
	}	
}

absorber::absorber(
	const boost::property_tree::ptree::value_type& absorbers_tree
):rectangle(absorbers_tree){}


