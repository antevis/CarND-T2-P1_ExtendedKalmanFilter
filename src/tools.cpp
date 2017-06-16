#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size() || estimations.size() == 0) {
        return rmse;
    }
    
    //accumulate squared diffs
    for(int i=0; i < estimations.size(); ++i){
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
    MatrixXd Hj(3,4);
    
    cout<<"x_state size: "<<x_state.size()<<"\n";
    
    for (int i = 0; i < x_state.size(); ++i) {
        
        cout<<"x_state("<<i<<"): "<< x_state(i)<<"\n";
    }
    
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //check division by zero
    if (px == 0 && py == 0) {
        return Hj;
    }
    
    //compute the Jacobian matrix
    const float pyxSqSum = pow(px, 2) + pow(py, 2);
    
    const float r0c0 = px / sqrt(pyxSqSum);
    const float r0c1 = py / sqrt(pyxSqSum);
    const float r0c2 = 0;
    const float r0c3 = 0;
    
    const float r1c0 = -py / pyxSqSum;
    const float r1c1 = px / pyxSqSum;
    const float r1c2 = 0;
    const float r1c3 = 0;
    
    const float r2c0 = py * (vx*py - vy*px)/pow(pyxSqSum, 2./3.);
    const float r2c1 = px * (vy*px - vx*py)/pow(pyxSqSum, 2./3.);
    const float r2c2 = px / sqrt(pyxSqSum);
    const float r2c3 = py / sqrt(pyxSqSum);
    
    Hj <<
        r0c0, r0c1, r0c2, r0c3,
        r1c0, r1c1, r1c2, r1c3,
        r2c0, r2c1, r2c2, r2c3;
    
    return Hj;
}

VectorXd Tools::cartesianToPolar(const VectorXd& xState){
    
//    const double THRESH = 0.0001;
    VectorXd polar(3);
    
    const double px = xState(0);
    const double py = xState(1);
    const double vx = xState(2);
    const double vy = xState(3);
    
    const double rho = sqrt(px * px + py * py);
    const double theta = atan2(py, px);
    const double rhoDot = ( px * vx + py * vy ) / rho;
    
    polar << rho, theta, rhoDot;
    return polar;
}

VectorXd Tools::polarToCartesian(const VectorXd& xState){
    
    VectorXd cartesian(4);
    
    const double rho = xState(0);
    const double phi = xState(1);
    const double drho = xState(2);
    
    const double px = rho * cos(phi);
    const double py = rho * sin(phi);
    const double vx = drho * cos(phi);
    const double vy = drho * sin(phi);
    
    cartesian << px, py, vx, vy;
    return cartesian;
}
