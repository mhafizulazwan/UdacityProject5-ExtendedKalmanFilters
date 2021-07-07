#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing covariance matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   * Done
   */
  
  //measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // measurement matrix - radar
  // will be redefined later
  Hj_ << 0, 0, 0, 0,
         0, 0, 0, 0,
         0, 0, 0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */
    
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    // create the state covariance matrix
    ekf_.P_ = MatrixXd(4, 4);;
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;
    
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;
    
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << 0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0,
               0, 0, 0, 0;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state. (Done)
      
      float ro = measurement_pack.raw_measurements_[0];  // get rho measurement
      float theta = measurement_pack.raw_measurements_[1]; // get theta measurement
      //float ro_dot = measurement_pack.raw_measurements_[2]; // get ro_dot measurement
        
      // convert to cartesian coordinates
      float px = cos(theta)*ro;
      float py = sin(theta)*ro;
      float vx = 0; //is it okay to use cos(theta)*ro_dot?
      float vy = 0; //is it okay to use sin(theta)*ro_dot?
      
      // initialize the state 
      ekf_.x_ << px,py,vx,vy;
      
      // all states can be measured and thus all covariance values are set to be lower (can be tuned)
      ekf_.P_(2,2) = 10; //vx
      ekf_.P_(3,3) = 10; //vy

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state. (Done)
      float px = measurement_pack.raw_measurements_[0]; // get Px position
      float py = measurement_pack.raw_measurements_[1]; // get Py position
      float vx = 0; // the velocity info is still unknown
      float vy = 0; // the velocity info is still unknown
      
      // initialize the state 
      ekf_.x_ << px,py,vx,vy;
      
      // speed can't be measured and thus the covariance values for speed are set to be higher (can be tuned)
      ekf_.P_(2,2) = 1000; //vx
      ekf_.P_(3,3) = 1000; //vy
      
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;   
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   * Done
   */
  
  //compute change of time (microseconds) and convert to seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0; 
  previous_timestamp_ = measurement_pack.timestamp_;
  
  // update the state transition matrix
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  
  // update the process noise covariance matrix
  float dt2 = dt*dt;
  float dt3 = dt*dt2;
  float dt4 = dt*dt3;
  int noise_ax = 9;
  int noise_ay = 9;
  
  ekf_.Q_(0,0) = (dt4/4)*noise_ax;
  ekf_.Q_(0,2) = (dt3/2)*noise_ax;
  ekf_.Q_(1,1) = (dt4/4)*noise_ay;
  ekf_.Q_(1,3) = (dt3/2)*noise_ay;
  ekf_.Q_(2,0) = (dt3/2)*noise_ax;
  ekf_.Q_(2,2) = dt2*noise_ax;
  ekf_.Q_(3,1) = (dt3/2)*noise_ay;
  ekf_.Q_(3,3) = dt2*noise_ay;
  
  
  ekf_.Predict();
  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   * Done
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates (Done)
    ekf_.R_ = MatrixXd(3,3);
    ekf_.H_ = MatrixXd(3,4);
    ekf_.R_ = R_radar_;
    Hj_ = tools.CalculateJacobian(ekf_.x_); //use predicted state to compute Jacobian matrix
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates (Done)
    ekf_.R_ = MatrixXd(2,2);
    ekf_.H_ = MatrixXd(2,4);
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
