#include "FusionEKF.h"
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

    noise_ax = 16.5;
    noise_ay = 11;

    H_laser_ << 1.0, 0, 0, 0,
                0, 1.0, 0, 0;

    float s = 0.0254;
    R_laser_ << s, 0,
                0, s;

    s = 0.1;
    R_radar_ << s, 0, 0,
                0, s, 0,
                0, 0, s;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    if (!is_initialized_) {
        Initialize(measurement_pack);
        is_initialized_ = true;
        return;
    }

    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    if (!IsValidMeasurement(measurement_pack))
        return;

    Predict(dt);

    Update(measurement_pack);
}

bool FusionEKF::IsValidMeasurement(const MeasurementPackage &measurement_pack) const {
    return !(measurement_pack.raw_measurements_[0] == 0 && measurement_pack.raw_measurements_[1] == 0);
}

void FusionEKF::Update(const MeasurementPackage &measurement_pack) {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        auto Hj = tools.CalculateJacobian(ekf_.x_);
        ekf_.UpdateRadar(measurement_pack.raw_measurements_, Hj, R_radar_);
    } else {
        // Laser updates
        ekf_.UpdateLazer(measurement_pack.raw_measurements_, H_laser_, R_laser_);
    }
}

void FusionEKF::Predict(double dt) {
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    ekf_.Q_ = tools.CalculateQ(dt, noise_ax, noise_ay);

    ekf_.Predict();
}

void FusionEKF::Initialize(const MeasurementPackage &measurement_pack) {
    // first measurement
    previous_timestamp_ = measurement_pack.timestamp_;

    MatrixXd F(4, 4);
    F <<    1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    MatrixXd P(4, 4);
    P <<    1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    MatrixXd Q(4, 4);

    VectorXd x(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        x = InitRadar(measurement_pack);

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        x = InitLazer(measurement_pack);
    }

    ekf_.Init(x, P, F, Q);
}

VectorXd FusionEKF::InitLazer(const MeasurementPackage &measurement_pack) {
    VectorXd x(4);
    x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    return x;
}

VectorXd FusionEKF::InitRadar(const MeasurementPackage &measurement_pack) {
    VectorXd x(4);
    float rho = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float rho_dot = measurement_pack.raw_measurements_[2];

    auto p = Tools::ToCartesian(rho, phi);
    auto v = Tools::ToCartesian(rho_dot, phi);
    x << p[0], p[1], v[0], v[1];
    return x;
}

