#include "kalman_filter.h"

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::UpdateLazer(const VectorXd &z, const MatrixXd &H, const MatrixXd &R) {
    VectorXd y = z - H * x_;
    Update(y, H, R);
}

void KalmanFilter::Update(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
    MatrixXd S = H * P_ * H.transpose() + R;
    MatrixXd K = P_ * H.transpose() * S.inverse();

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H) * P_;
}

void KalmanFilter::UpdateRadar(const VectorXd &z, const MatrixXd &H, const MatrixXd &R) {
    double px = x_[0];
    double py = x_[1];
    double vx = x_[2];
    double vy = x_[3];

    float rho = sqrt(pow(px, 2) + pow(py, 2));
    float phi = atan(py/px);
    float rho_dot = (px*vx + py*vy) / rho;

    VectorXd pred(3,1);
    pred << rho, phi, rho_dot;

    VectorXd y = z - pred;
    Update(y, H, R);
}
