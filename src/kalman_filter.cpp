#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {
  
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  int size = x_.size();
  I_ = MatrixXd::Identity(size, size);

}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */

  // y = current measurement - predicted measurement
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();

  //new state info
  x_ = x_ + K * y;
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  float rho = sqrt(px * px + py * py);
  float phi = atan2(py, px);
  float rho_dot = (px * vx + py * vy) / rho;

  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rho_dot;

  VectorXd y = z - z_pred;

  // One other important point when calculating y with radar sensor data: 
  // the second value in the polar coordinate vector is the angle ϕ. 
  // You'll need to make sure to normalize ϕ in the y vector so that 
  // its angle is between −pi and pi; in other words, add or subtract 2pi 
  // from ϕ until it is between −pi and pi.
  while (y[1] < -M_PI) {
    y[1] += 2 * M_PI;
  }
  while (y[1] > M_PI) {
    y[1] -= 2 * M_PI;
  }

  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();

  x_ = x_ + K * y;
  P_ = (I_ - K * H_) * P_;
}
