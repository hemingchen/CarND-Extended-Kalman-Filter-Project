#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // Calculate predicted radar measurements
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  float rho;
  float phi;
  float rho_dot;

  // Ensure px and py are big enough therefore no division by zero in calculating Hj
  float p_threshold = 0.0001;
  if(fabs(px) < p_threshold){
    cout << "px too small = " << px << endl;
    px = p_threshold;
    cout << "px bumped to = " << px << endl;
  }

  if(fabs(py) < p_threshold){
    cout << "py too small = " << py << endl;
    py = p_threshold;
    cout << "py bumped to = " << py << endl;
  }

  rho = sqrt(px * px + py * py);
  phi = atan2(py, px); // atan of y/x has to be in [-pi,+pi].
  rho_dot = (px * vx + py * vy) / rho;

  // Get predicted measurements
  VectorXd hx(3);
  hx << rho, phi, rho_dot;

  // Calculate measurement update
  VectorXd y = z - hx;

  // Ensure phi is eventually normalized into [-pi, +pi]
  float PI = 3.1415926;
  cout << "before regulation phi in y  = " << phi / PI << " pi" << endl;
  while (fabs(y(1)) > PI) {
    if (y(1) > 0) {
      y(1) -= 2 * PI;
    } else {
      y(1) += 2 * PI;
    }
  }
  cout << "after regulation  phi in y  = " << phi / PI << " pi" << endl;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
