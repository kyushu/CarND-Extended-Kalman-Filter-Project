#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
	is_initialized_ = false;

	previous_timestamp_ = 0;

	// initializing matrices
	R_laser_ = MatrixXd(2, 2);
	R_radar_ = MatrixXd(3, 3);
	H_laser_ = MatrixXd(2, 4);
	Hj_ = MatrixXd(3, 4);

	//measurement covariance matrix - laser
	R_laser_ << 0.0225, 0,
				0, 0.0225;

	//measurement covariance matrix - radar
	R_radar_ << 0.09, 0,      0,
				0,    0.0009, 0,
				0,    0,      0.09;

	/**
	 TODO:
	* Finish initializing the FusionEKF.
	* Set the process and measurement noises
	*/
	noise_ax = 9.0f;
	noise_ay = 9.0f;

	H_laser_ << 1, 0, 0, 0,
			  	0, 1, 0, 0;

	ekf_.F_ = MatrixXd(4,4);
	ekf_.F_ << 1, 0, 1, 0,
				0, 1, 0, 1,
				0, 0, 1, 0,
				0, 0, 0, 1;

	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

	ekf_.Q_ = MatrixXd(4,4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


	/*****************************************************************************
	 *  Initialization
	 ****************************************************************************/
	if (!is_initialized_) {
		/**
		TODO:
			* Initialize the state ekf_.x_ with the first measurement.
			* Create the covariance matrix.
			* Remember: you'll need to convert radar from polar to cartesian coordinates.
		*/

		// first measurement
		cout << "EKF: " << endl;
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << 1, 1, 1, 1;

		// measurement(x) = [initial_x, initial_y,velocity_x, velocity_y]

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			 Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = measurement_pack.raw_measurements_[0];
			float phi = measurement_pack.raw_measurements_[1];
			float x = rho * cosf(phi);
			float y = rho * sinf(phi);
			
			// i don't know hot to calculate velocity_x and velocity_y according to
			// rho, phi, rho_dot, so i set vx = vy = 0
			ekf_.x_ << x, y, 0, 0;
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			 Initialize state.
			*/
			//set the state with the initial location and zero velocity
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
			previous_timestamp_ = measurement_pack.timestamp_;
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}
	

	/*****************************************************************************
	 *  Prediction
	 ****************************************************************************/

	/**
	 TODO:
		* Update the state transition matrix F according to the new elapsed time.
		- Time is measured in seconds.
		* Update the process noise covariance matrix.
		* Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
	*/
	// Calculate the new elapsed time
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	// update previous_timestamp 
	previous_timestamp_ = measurement_pack.timestamp_;

	// Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	//2. Set the process covariance matrix Q
	// Given noise_ax and noise_ay
	double dt2 = dt*dt;
	double dt3 = dt2*dt;
	double dt4 = dt3*dt;

	double vx1 = dt4/4.0 * noise_ax;
	double vx2 = dt3/2.0 * noise_ax;
	double vx3 = dt2*noise_ax;

	double vy1 = dt4/4.0 * noise_ay;
	double vy2 = dt3/2.0 * noise_ay;
	double vy3 = dt2*noise_ay;

	ekf_.Q_ << vx1,   0, vx2,   0,
			     0, vy1,   0, vy2,
			   vx2,   0, vx3,   0,
			     0, vy2,   0, vy3;                 


	// ekf_.Q_(0, 0) = vx1;
	// ekf_.Q_(2, 2) = vx3;
	// ekf_.Q_(0, 2) = vx2;
	// ekf_.Q_(2, 0) = vx2;

	// ekf_.Q_(1, 1) = vy1;
	// ekf_.Q_(3, 3) = vy3;
	// ekf_.Q_(1, 3) = vy2;
	// ekf_.Q_(3, 1) = vy2;

	ekf_.Predict();
	/*****************************************************************************
	 *  Update
	 ****************************************************************************/

	/**
	 TODO:
		* Use the sensor type to perform the update step.
		* Update the state and covariance matrices.
	*/

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		ekf_.R_ = R_radar_;
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.H_ = Hj_;
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	} else {
		// Laser updates
		ekf_.R_ = R_laser_;
		ekf_.H_ = H_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
