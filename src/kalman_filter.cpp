#include "kalman_filter.h"
#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {
    
    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    
    
    //measurement covariance matrix - laser
    R_laser_ <<
        0.0225, 0,
        0, 0.0225;
    
    //measurement covariance matrix - radar
    R_radar_ <<
        0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
    
    H_laser_ <<
        1, 0, 0, 0,
        0, 1, 0, 0;
    
    
    F_ = Eigen::MatrixXd(4, 4);
    F_ <<
        1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
    
    P_ = Eigen::MatrixXd(4, 4);
    P_ <<
        1, 0, 0,    0,
        0, 1, 0,    0,
        0, 0, 1000, 0,
        0, 0, 0,    1000;
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
    
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

Eigen::MatrixXd KalmanFilter::KMatrix() {
    MatrixXd h_t = H_.transpose(); // used more than once, so avoiding repeated calculations
    MatrixXd s = H_ * P_ * h_t + R_;
    MatrixXd K = P_ * h_t * s.inverse();
    
    return K;
}

Eigen::MatrixXd KalmanFilter::QMatrix(float dt, float noise_ax, float noise_ay) {
    
    Eigen::MatrixXd q = Eigen::MatrixXd(4, 4);
    
    const float r0c0 = pow(dt, 4) / 4 * noise_ax;
    const float r0c1 = 0;
    const float r0c2 = pow(dt,3) / 2 * noise_ax;
    const float r0c3 = 0;
    
    const float r1c0 = 0;
    const float r1c1 = pow(dt, 4) / 4 * noise_ay;
    const float r1c2 = 0;
    const float r1c3 = pow(dt, 3) / 2 * noise_ay;
    
    const float r2c0 = r0c2;
    const float r2c1 = 0;
    const float r2c2 = pow(dt, 2) * noise_ax;
    const float r2c3 = 0;
    
    const float r3c0 = 0;
    const float r3c1 = r1c3;
    const float r3c2 = 0;
    const float r3c3 = pow(dt, 2) * noise_ay;
    
    q <<
        r0c0, r0c1, r0c2, r0c3,
        r1c0, r1c1, r1c2, r1c3,
        r2c0, r2c1, r2c2, r2c3,
        r3c0, r3c1, r3c2, r3c3;
    
    return q;
}

void KalmanFilter::SetQMatrix(float dt, float noise_ax, float noise_ay) {
    Q_ = QMatrix(dt, noise_ax, noise_ay);
}

void KalmanFilter::SetFMatrix(float dt) {
    F_(0,2) = dt;
    F_(1,3) = dt;
}

void KalmanFilter::Update(const MeasurementPackage& package) {
    
    // Update
    if (package.sensor_type_ == MeasurementPackage::RADAR) {
        
        Eigen::VectorXd cartesianMeas = aux_.PolarToCartesian(package.raw_measurements_);
        
        H_ = aux_.CalculateJacobian(cartesianMeas);
        R_ = R_radar_;
        UpdateWithRadar(package.raw_measurements_);
        
        
    } else {
        // Laser updates
        H_ = H_laser_;
        R_ = R_laser_;
        UpdateWithLidar(package.raw_measurements_);
        
    }
}

void KalmanFilter::ComputeNewEstimate(const Eigen::MatrixXd& K, const Eigen::VectorXd& y) {
    
    x_ = x_ + (K * y);
    const long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateWithLidar(const VectorXd& z) {
    
    VectorXd y = z - H_ * x_;
    
    //new estimate
    ComputeNewEstimate(KMatrix(), y);
}

void KalmanFilter::UpdateWithRadar(const VectorXd& z) {
    
    VectorXd pred = aux_.CartesianToPolar(x_);
    VectorXd y = z - pred;
    
    //Normalizing phi as suggested by Slack community:
    //https://carnd.slack.com/files/udasuburb/F5R5S9XEU/pasted_image_at_2017_06_10_02_39_pm.png
    y(1) = fmod(y(1), M_PI);
    
    //new estimate
    ComputeNewEstimate(KMatrix(), y);
}
