#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

const float EPS = 0.0001;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
															const vector<VectorXd> &ground_truth) {

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];
		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj = MatrixXd::Zero(3, 4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// printf("px: %f\n", px);
	// printf("py: %f\n", py);
	// printf("vx: %f\n", vx);
	// printf("vy: %f\n", vy);
	
	//compute the Jacobian matrix
	float c1 = px*px+py*py;
	
	//check division by zero
	if(c1 < EPS) {
		return Hj;
	}
	
	float c2 = sqrt(c1);
	float c3 = (c1*c2);
		
	float ro_dpx = px / c2;
	float ro_dpy = py / c2;
	float ro_dvx = 0;
	float ro_dvy = 0;

	float theta_dpx = -py / c1;
	float theta_dpy = px / c1;
	float theta_dvx = 0;
	float theta_dvy = 0;

	printf("%f\n", vy*px - vx*py);

	float rop_dpx = py*(vx*py - vy*px) / c3;
	float rop_dpy = px*(vy*px - vx*py) / c3;
	float rop_dvx = px / c2;
	float rop_dvy = py / c2;

	Hj << ro_dpx,    ro_dpy,    ro_dvx,    ro_dvy,
		  theta_dpx, theta_dpy, theta_dvx, theta_dvy,
		  rop_dpx,   rop_dpy,   rop_dvx,   rop_dvy;

	return Hj;
}
