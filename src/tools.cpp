#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd>& estimations,
                              const vector<VectorXd>& ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size() || estimations.size() == 0) {
        return rmse;
    }
    
    //accumulate squared diffs
    for(int i=0; i < estimations.size(); ++i) {
        VectorXd currDiff = estimations[i] - ground_truth[i];
        
        currDiff = currDiff.array()*currDiff.array();
        
        rmse += currDiff;
    }
    
    //calculate the mean
    rmse /= estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd h_j(3,4);
    
    //recover state parameters
    const float px = x_state(0);
    const float py = x_state(1);
    const float vx = x_state(2);
    const float vy = x_state(3);
    
    // check for division by zero
    // denominator becomes zero iif both px and py are zero.
    if (px == 0 && py == 0) {
        return h_j;
    }
    
    //compute the Jacobian matrix
    const float py_px_sq_sum = pow(px, 2) + pow(py, 2);
    
    const float r0c0 = px / sqrt(py_px_sq_sum);
    const float r0c1 = py / sqrt(py_px_sq_sum);
    const float r0c2 = 0;
    const float r0c3 = 0;
    
    const float r1c0 = -py / py_px_sq_sum;
    const float r1c1 = px / py_px_sq_sum;
    const float r1c2 = 0;
    const float r1c3 = 0;
    
    const float r2c0 = py * (vx*py - vy*px)/pow(py_px_sq_sum, 2./3.);
    const float r2c1 = px * (vy*px - vx*py)/pow(py_px_sq_sum, 2./3.);
    const float r2c2 = px / sqrt(py_px_sq_sum);
    const float r2c3 = py / sqrt(py_px_sq_sum);
    
    h_j <<
        r0c0, r0c1, r0c2, r0c3,
        r1c0, r1c1, r1c2, r1c3,
        r2c0, r2c1, r2c2, r2c3;
    
    return h_j;
}

//Credits to Mithi Sevilla: https://github.com/mithi/fusion-ekf/blob/master/src/tools.cpp
VectorXd Tools::CartesianToPolar(const VectorXd& measurement) {
    
    VectorXd polar(3);
    
    const double px = measurement(0);
    const double py = measurement(1);
    const double vx = measurement(2);
    const double vy = measurement(3);
    
    const double rho = sqrt(px * px + py * py);
    const double phi = atan2(py, px);
    const double rho_dot = ( px * vx + py * vy ) / rho;
    
    polar << rho, phi, rho_dot;
    return polar;
}

//Credits to Mithi Sevilla: https://github.com/mithi/fusion-ekf/blob/master/src/tools.cpp
VectorXd Tools::PolarToCartesian(const VectorXd& measurement) {
    
    VectorXd cartesian(4);
    
    const double rho = measurement(0);
    const double phi = measurement(1);
    const double rho_dot = measurement(2);
    
    const double px = rho * cos(phi);
    const double py = rho * sin(phi);
    const double vx = rho_dot * cos(phi);
    const double vy = rho_dot * sin(phi);
    
    cartesian << px, py, vx, vy;
    return cartesian;
}
