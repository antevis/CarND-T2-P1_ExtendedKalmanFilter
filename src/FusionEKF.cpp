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
}

// Destructor.
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage& measurement_pack) {

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
            Eigen::VectorXd init_state = tools.PolarToCartesian(measurement_pack.raw_measurements_);
            
            ekf_.x_(0) = init_state(0);
            ekf_.x_(1) = init_state(1);
            
        } else {
            
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
    
    //1. Updating the F matrix so that the time is integrated
    ekf_.SetFMatrix(dt);
    
    //2. Updating the Q matrix
    ekf_.SetQMatrix(dt, noise_ax_, noise_ay_);
    
    //3.
    ekf_.Predict();

    //4.
    ekf_.Update(measurement_pack);

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
