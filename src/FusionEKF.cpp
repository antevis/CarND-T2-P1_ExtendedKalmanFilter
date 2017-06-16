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
    R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

    
    H_laser_ <<
        1, 0, 0, 0,
        0, 1, 0, 0;
    
    ekf_.F_ = Eigen::MatrixXd(4, 4);
    ekf_.F_ <<
        1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;
    
    ekf_.P_ = Eigen::MatrixXd(4, 4);
    ekf_.P_ <<
        1, 0, 0,    0,
        0, 1, 0,    0,
        0, 0, 1000, 0,
        0, 0, 0,    1000;
}

// Destructor.
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {
        
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            
            //Convert radar from polar to cartesian coordinates and initialize state.
            Eigen::VectorXd init_state = aux_.PolarToCartesian(measurement_pack.raw_measurements_);
            
            ekf_.x_(0) = init_state(0);
            ekf_.x_(1) = init_state(1);
            
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            
            //Initialize state.
            ekf_.x_(0) = measurement_pack.raw_measurements_[0];
            ekf_.x_(1) = measurement_pack.raw_measurements_[1];
        }
        
        previous_timestamp_ = measurement_pack.timestamp_;
      
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/
    //compute the time elapsed between the current and previous measurements
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // TODO: YOUR CODE HERE
    //1. Modify the F matrix so that the time is integrated
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;
    
    //2. Set the process covariance matrix Q
    ekf_.Q_ = Eigen::MatrixXd(4, 4);
    const float r0c0 = pow(dt, 4) / 4 * noise_ax_;
    const float r0c1 = 0;
    const float r0c2 = pow(dt,3) / 2 * noise_ax_;
    const float r0c3 = 0;
    
    const float r1c0 = 0;
    const float r1c1 = pow(dt, 4) / 4 * noise_ay_;
    const float r1c2 = 0;
    const float r1c3 = pow(dt, 3) / 2 * noise_ay_;
    
    const float r2c0 = r0c2;
    const float r2c1 = 0;
    const float r2c2 = pow(dt, 2) * noise_ax_;
    const float r2c3 = 0;
    
    const float r3c0 = 0;
    const float r3c1 = r1c3;
    const float r3c2 = 0;
    const float r3c3 = pow(dt, 2) * noise_ay_;
    
    ekf_.Q_ <<
        r0c0, r0c1, r0c2, r0c3,
        r1c0, r1c1, r1c2, r1c3,
        r2c0, r2c1, r2c2, r2c3,
        r3c0, r3c1, r3c2, r3c3;

    ekf_.Predict();

    // Update
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            
        Eigen::VectorXd cartesianMeas = aux_.PolarToCartesian(measurement_pack.raw_measurements_);
        
        ekf_.H_ = aux_.CalculateJacobian(cartesianMeas);
        ekf_.R_ = R_radar_;
        ekf_.UpdateWithRadar(measurement_pack.raw_measurements_);
        
        
    } else {
    // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.UpdateWithLidar(measurement_pack.raw_measurements_);
        
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
