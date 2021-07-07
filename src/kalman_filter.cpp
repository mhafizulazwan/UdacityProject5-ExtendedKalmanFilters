#include "kalman_filter.h"
#include <cmath>
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   * TODO: predict the state
   * Done
   */
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_; 
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   * Done
   */
  VectorXd y = z - H_*x_;
  MatrixXd S = H_*P_*H_.transpose() + R_;
  MatrixXd K = P_*H_.transpose()*S.inverse();
  
  x_ = x_ + (K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_)*P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   * Done
   */
  // get the predicted state
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  //calculate non-linear function h(x')
  VectorXd hx(3);
  hx(0) = sqrt(px*px + py*py); 
  hx(1) = atan2(py,px);
  hx(2) = (px*vx + py*vy)/(sqrt(px*px + py*py));
  
  VectorXd y = z - hx; //y becomes the error of rho, phi, and rho_dot
  // normalize y(1) between -PI and PI
  while (y(1) < -M_PI || y(1) > M_PI){
    if (y(1) < -M_PI){
      y(1) += 2*M_PI;
    } else {
      y(1) -= 2*M_PI;
    }
  }
  
  MatrixXd S = H_*P_*H_.transpose() + R_; // Hj (Jacobian matrix) is used instead of nonlinear function h(x)
  MatrixXd K = P_*H_.transpose()*S.inverse();
  
  x_ = x_ + (K*y); // the new state is 4x1 instead of 3x1
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_)*P_;
}
