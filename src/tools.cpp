#include <iostream>
#include <cmath>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
    TODO:
      * Calculate the RMSE here.
    */
    VectorXd rmse(4);
    rmse << 0.0 , 0.0, 0.0, 0.0;

    if (estimations.size() ==0 || estimations.size() != ground_truth.size())
        return rmse;

    for (int i =0; i < estimations.size(); ++i)
    {
        VectorXd diff = estimations[i] - ground_truth[i];
        diff = diff.array() * diff.array();
        rmse += diff;
    }

    rmse /= estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }

    //compute the Jacobian matrix
    Hj <<   (px/c2), (py/c2), 0, 0,
            -(py/c1), (px/c1), 0, 0,
            py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

    return Hj;
}

MatrixXd Tools::CalculateQ(const double dt, float noise_ax, float noise_ay) {
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    MatrixXd Q(4, 4);

    Q <<    dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
            0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
            dt_3/3*noise_ax,    0, dt_2*noise_ax, 0,
            0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

    return Q;
}

VectorXd Tools::ToCartesian(float rho, float phi) {
    VectorXd result(2,1);
    result << rho * cos(phi), rho * sin(phi);
    return result;
}
