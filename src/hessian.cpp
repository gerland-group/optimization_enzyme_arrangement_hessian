#include "hessian.h"

void hessian(
	Eigen::MatrixXd& H,
	Eigen::VectorXd& U,
	const mutant& before,
	const double alpha,
	const int Nsites,
	const std::vector<int>& pos
){
	Eigen::MatrixXd O = Eigen::MatrixXd(Nsites,Nsites);
	//cout << "alpha Mt e" << endl;
	//cout << (alpha*before[0].inv_M.transpose()*before[0].E) << endl;
	//cout << "1-alpha Mt e" << endl;
	//cout << U - (alpha*before[0].inv_M.transpose()*before[0].E) << endl;
		
	//cout << "transpose" << endl;
	//cout << ((U - (alpha*before[0].inv_M.transpose()*before[0].E)).transpose()) << endl;
	
	O = before.rho*((U + (alpha*before.inv_M.transpose()*before.E)).transpose());
	
	//cout << "O" << endl;
	//cout << O << endl;
	//cout << O(0,0) << endl;
	
	for(unsigned int j=0; j<pos.size(); j++){
		for(unsigned int k=0; k<pos.size(); k++){
			H(j,k) = alpha*alpha*(  before.inv_M(pos[j],pos[k])*O(pos[k],pos[j])  +  before.inv_M(pos[k],pos[j])*O(pos[j],pos[k]) );
		}
	}
	
	O.resize(0,0);
	
}


void solve_de(
	Eigen::VectorXd& de,
	Eigen::MatrixXd& H
){
	//Solve the system H*de=1 with LU decomposition
	Eigen::VectorXd id=Eigen::VectorXd::Ones(de.size());
	
	de=H.partialPivLu().solve(id);
	
	//Normalization
	double norm = de.sum();
	de/=norm;
	
	/*cout << "Vector de has components:" << endl;
	for(int i=0;i<de.size();i++){
		cout << de(i) << endl;
	}*/
}
