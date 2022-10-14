//
//  FEM.h
//  CompHydro
//
//  Created by Olga Diamanti on 27/02/15.
//  Copyright (c) 2015 Olga Diamanti. All rights reserved.
//

#ifndef __CompHydro__FEM__
#define __CompHydro__FEM__

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <vector>
#include <string>

using namespace std;
using namespace Eigen;

class FEM
{
public:
  FEM();
  FEM(const Eigen::MatrixXd &_Vf,
      const Eigen::MatrixXi &_Ff,
      const int _inner_iter,
      const Eigen::VectorXd &_young_modulus,
      const double _poisson_ratio,
      const bool _fem_symmetric = false);

  bool solve(const Eigen::VectorXi &b,
             Eigen::MatrixXd &Vf_uv, double current_step=1);
  Eigen::Matrix2d S;
  Eigen::Matrix2d I;

//protected:
  int inner_iter;
  double poisson_ratio;

  int matrix_type;
  double current_step;
  Eigen::MatrixXd Bm;
  Eigen::VectorXd W;
  const Eigen::MatrixXd &Vf;
  const Eigen::MatrixXi &Ff;
  const Eigen::VectorXd &young_modulus;
  Eigen::VectorXd mu, lambda;
  Eigen::MatrixXd FNf;

  void flatten_all(const Eigen::MatrixXd &Vf_uv,
                   Eigen::MatrixXd &Xflat,
                   Eigen::MatrixXd &xflat);

  void precompute (const Eigen::MatrixXd &X);

  double elastic_energy(const Eigen::MatrixXd &x);
  void elastic_forces(const Eigen::MatrixXd &x,
                      Eigen::VectorXd &fe,
                      bool do_derivative = false,
                      Eigen::SparseMatrix<double> &dfedx = *(Eigen::SparseMatrix<double>*)NULL,
                      double current_step=1);
  Eigen::Vector4d deformation_gradients(const Eigen::MatrixXd &x,
                                        const int fi,
                                        bool do_derivative = false,
                                        Eigen::Matrix<double,4,6> &dFdx = *(Eigen::Matrix<double,4,6>*)NULL);
  Eigen::Vector4d corotational_elasticity_stress_tensor(const Eigen::MatrixXd &x,
                                                        const int fi,
                                                        bool do_derivative,
                                                        Eigen::Matrix<double,4,6> &dPdx = *(Eigen::Matrix<double,4,6>*)NULL,
                                                        double current_step=1);
  Eigen::MatrixXd xflat;

};

#endif /* defined(__CompHydro__FEM__) */
