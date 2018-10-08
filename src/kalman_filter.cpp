#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;


// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.
const double PI  =3.141592653589793238463;
const float EPS = 0.0001;


KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
												MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict() {
	/**
	TODO:
		* predict the state
	*/
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	/**
	TODO:
		* update the state by using Kalman Filter equations
	*/
	/* 
	 * x_(prediction) = [px, py, vx, vy]
	 * H (Transition matrix) = [ 1,0,0,0
     *                           0,1,0,0 ] 
	 * z_pred = H*x_;
	 */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	TODO:
		* update the state by using Extended Kalman Filter equations
	*/
	/*
	 * x_(prediction) = [px, py, vx, vy]
	 * h(x_) (Transition function) = [ sqrt(px^2+py^),
	 *                                 arctan(py/px) ,
	 *                                 pxvx+pyvy / sqrt(px^2+py^) ]
	 * z_pred = h(x_)
	 */
	double px = x_[0];
	double py = x_[1];
	double vx = x_[2];
	double vy = x_[3];

	double rho = sqrt(px*px+py*py);
	double phi = 0;
	double rho_dot = 0;

	if(rho > EPS)
	{
		phi = atan2(py, px);
		rho_dot = (px*vx+py*vy) / rho;
	}

	VectorXd z_pred(3);
	z_pred << rho, phi, rho_dot;

	VectorXd y = z - z_pred;

	// Restrict the range of theta is between -PI and PI 
	while(y[1] < PI) {
		y[1] += (2*PI);
	}
	while(y[1] > PI) {
		y[1] -= (2*PI);
	}

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
