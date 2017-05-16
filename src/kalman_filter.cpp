#include "kalman_filter.h"
# define M_PI           3.14159265358979323846  /* pi */

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
    /*
     * KF Prediction step
     */
    x_ = F_ * x_;// + u;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /*
     * KF Measurement update step
     */
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    //new state
    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /*
     * KF Measurement update step
     */
     // Non linearity
    float r = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    float tet = atan2(x_(1), x_(0));
    float r_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / r;
    VectorXd h = VectorXd(3);
    h << r, tet, r_dot;

    VectorXd y = z - h;
    // we waork with the angle in the range of [-pi,pi)
    y(1) = fmod(y(1)+M_PI, (2.0*M_PI)) - M_PI;

    // same as in Update
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);

    //new state
    x_ = x_ + (K * y);
    P_ = (I - K * H_) * P_;
}
