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

    /**
    * Updates the state by using standard Kalman Filter equations
    * @param z The measurement at k+1
    */
    void UpdateWithLidar(const Eigen::VectorXd &z);

    /**
    * Updates the state by using Extended Kalman Filter equations
    * @param z The measurement at k+1
    */
    void UpdateWithRadar(const Eigen::VectorXd &z);
    
    void ComputeNewEstimate(const Eigen::MatrixXd& K, const Eigen::VectorXd& y);
    
    Eigen::MatrixXd KMatrix();
    Eigen::MatrixXd QMatrix(float dt, float noise_ax, float noise_ay);
    void SetQMatrix(float dt, float noise_ax, float noise_ay);
    void SetFMatrix(float dt);
    
    void Update(const MeasurementPackage& package);
    
    //To conveniently invoke helper functions
    Tools aux_;

};

#endif /* KALMAN_FILTER_H_ */
