#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

    /**
    * A helper method to calculate Jacobians.
    */
    MatrixXd CalculateJacobian(const VectorXd &x_state);

    /**
    * A helper function to calculate covariance noise
    */
    MatrixXd CalculateQ(const double dt, float noise_ax, float noise_ay);

    static VectorXd ToCartesian(float rho, float phi);

};

#endif /* TOOLS_H_ */
