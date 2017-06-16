#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"
#include "measurement_package.h"
#include "tools.h"

class KalmanFilter {
public:

    // state vector
    Eigen::VectorXd x_;

    // state covariance matrix
    Eigen::MatrixXd P_;

    // state transition matrix
    Eigen::MatrixXd F_;

    // process covariance matrix
    Eigen::MatrixXd Q_;

    // measurement matrix
    Eigen::MatrixXd H_;

    // measurement covariance matrix
    Eigen::MatrixXd R_;
    
    //Sensors data
    Eigen::MatrixXd R_laser_;
    Eigen::MatrixXd R_radar_;
    Eigen::MatrixXd H_laser_;

    /**
    * Constructor
    */
    KalmanFilter();

    /**
    * Destructor
    */
    virtual ~KalmanFilter();

    /**
    * Prediction Predicts the state and the state covariance
    * using the process model
    * @param delta_T Time between k and k+1 in s
    */
    void Predict();
    
    // Part of the Update phase. Computes x_ and P_
    void ComputeNewEstimate(const Eigen::MatrixXd& K, const Eigen::VectorXd& y);
    
    // Computes K matrix
    Eigen::MatrixXd KMatrix();
    // Updates Q matrix
    Eigen::MatrixXd QMatrix(float dt, float noise_ax, float noise_ay);
    void SetQMatrix(float dt, float noise_ax, float noise_ay);
    
    // Updates F matrix
    void SetFMatrix(float dt);
    
    /**
     * Updates the state by using Kalman Filter equations
     * @param package The measurement package at k+1
     */
    void Update(const MeasurementPackage& package);
    
    // Computes y
    Eigen::VectorXd RadarLoss(const Eigen::VectorXd& z);
    Eigen::VectorXd LidarLoss(const Eigen::VectorXd& z);
    
    //To conveniently invoke helper functions
    Tools aux_;

};

#endif /* KALMAN_FILTER_H_ */
