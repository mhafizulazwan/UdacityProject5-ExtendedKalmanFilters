#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   * Done
   */
  
  //initialize RMSE
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  //accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i){
    
    VectorXd c1 = estimations[i] - ground_truth[i];
    VectorXd c2 = c1.array()*c1.array();
  	rmse += c2;
  }
  
  //calculate the mean
  rmse = rmse/estimations.size();
  
  //calculate the square root
  rmse = rmse.array().sqrt();
  
  //return the results
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   * Done
   */
  
  //initialize Jacobian matrix
  MatrixXd Hj(3,4);
  
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  
  float c1 = px*px + py*py;
  float c2 = sqrt(c1);
  float c3 = sqrt(c1*c1*c1);
  
  //check division by zero
  if (fabs(c1) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }
  
  //compute Jacobian matrix
  Hj(0,0) = px/c2; Hj(0,1) = py/c2; Hj(0,2) = 0; Hj(0,3) = 0;
  Hj(1,0) = -py/c1; Hj(1,1) = px/c1; Hj(1,2) = 0; Hj(1,3) = 0;
  Hj(2,0) = py*(vx*py - vy*px)/c3; Hj(2,1) = px*(vy*px - vx*py)/c3; Hj(2,2) = px/c2; Hj(2,3) = py/c2;
  
  return Hj;
}
