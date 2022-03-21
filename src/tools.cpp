#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd result(4);
  result << 0, 0, 0, 0;
  
  for (unsigned int i=0; i<estimations.size(); ++i) {
    VectorXd temp = estimations[i]-ground_truth[i];
    temp = temp.array()*temp.array();
    result += temp;
  }
  
  result = result/estimations.size();
  result = result.array().sqrt();
  
  return result;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd J(3,4);
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  if (fabs(px) < 0.001 and fabs(py) < 0.001) {
    px = (px > 0) ? 0.001 : -0.001;
    py = (py > 0) ? 0.001 : -0.001;
  }

  float c1 = px*px+py*py;
  c1 = (c1 < 0.0000001) ? 0.0000001 : c1;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  J << (px/c2),              (py/c2),                0,    0,
        -(py/c1),             (px/c1),                0,     0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return J;
}
