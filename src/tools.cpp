#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
   // Calculate the RMSE
	VectorXd rmse(4);
		rmse << 0., 0., 0., 0.;

		if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
			std::cout << "Invalid estimation or ground_truth data" << std::endl;
			return rmse;
		}

		for (unsigned int i = 0; i < estimations.size(); i++) {
			VectorXd residual = estimations[i] - ground_truth[i];
			residual = residual.array()*residual.array();
			rmse += residual;
		}

		rmse = rmse / estimations.size();
		rmse = rmse.array().sqrt();
		return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    // Calculate a Jacobian
	MatrixXd Hj(3, 4);
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	
	//pre-compute terms
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check div/0
	if (fabs(c1) < 0.00001) {
		std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
		return Hj;
	}

	//compute Jacobian
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;

	return Hj;

  
}
