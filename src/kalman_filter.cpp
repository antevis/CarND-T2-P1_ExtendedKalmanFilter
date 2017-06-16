#include "kalman_filter.h"

#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
    
    x_ = F_ * x_;
//    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    
    VectorXd y = z - H_ * x_;
    MatrixXd hT = H_.transpose();
    MatrixXd s = H_ * P_ * hT + R_;
    MatrixXd K = P_ * hT * s.inverse();
    
    //new estimate
    x_ = x_ + (K * y);
    const long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
    Tools t;
    
    VectorXd pred = t.cartesianToPolar(x_);
    
    VectorXd y = z - pred;
    
    y(1) = fmod(y(1), M_PI);
    
    MatrixXd hT = H_.transpose();
    MatrixXd s = H_ * P_ * hT + R_;
    MatrixXd K = P_ * hT * s.inverse();
    
    //new estimate
    x_ = x_ + (K * y);
    const long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
