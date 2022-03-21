#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   * TODO: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * (F_.transpose()) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd h = H_ * x_;
  VectorXd dh = z - h;
  UpdateWithdh(dh);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  double rho     = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double theta   = atan2(x_(1), x_(0));
  double rho_v   = (fabs(rho) < 0.001) ? 0 : (x_(0)*x_(2) + x_(1)*x_(3)) / rho;

  VectorXd h = VectorXd(3);
  h << rho, theta, rho_v;
  VectorXd dh = z - h;
  UpdateWithdh(dh);
}

void KalmanFilter::UpdateWithdh(VectorXd &y) {
  
  MatrixXd Ht  = H_.transpose();
  MatrixXd S   = H_ * P_ * Ht + R_;
  MatrixXd K   = P_ * Ht * (S.inverse());

//   if (x_[0] < -4.9 and x_[0] > -5.4) {
//     cout << "$$$$$$$$$$$$$$" << endl;
//     cout << "S: " << S << endl;
//     cout << "K: " << K << endl;
//     cout << "y: " << y << endl;
//     cout << "$$$$$$$$$$$$$$" <<endl;
//   }
  
  while (y[1] > M_PI)  y[1] -= 2.0 * M_PI;
  while (y[1] < -M_PI) y[1] += 2.0 * M_PI;
  
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
