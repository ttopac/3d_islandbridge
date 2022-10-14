
//
//  FEM.cpp
//  CompHydro
//
//  Created by Olga Diamanti on 27/02/15.
//  Copyright (c) 2015 Olga Diamanti. All rights reserved.
//

#include "FEM.h"
#include "Dipping.h"

#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/procrustes.h>
#include <igl/polar_svd.h>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/rotation_matrix_from_directions.h>
#include <igl/per_face_normals.h>
#include <vector>
#include <fstream>

#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

using namespace Eigen;
using namespace std;

#include <iostream>
#include <fstream>
#include <igl/matlab_format.h>
void dRdF_func(const Matrix2d &R,
               const Matrix2d &S,
               Matrix4d &dRdF)
{
  for (int i=0; i<2; ++i)
  {
    Vector2d ei; ei.setZero(); ei(i) = 1;
    for (int j=0; j<2; ++j)
    {
      Vector2d ej; ej.setZero(); ej(j) = 1;
      Matrix2d M = R.transpose()*ei*ej.transpose();
      M = M-M.transpose().eval();
      double dbl_skew = M(0,1);
      double x = dbl_skew / S.trace();
      Matrix2d dRdFij;
      dRdFij.row(0) = x*R.row(1);
      dRdFij.row(1) = -x*R.row(0);

      int col = i*2+j;
      dRdF.col(col) << dRdFij.row(0).transpose(), dRdFij.row(1).transpose();
    }
  }
}


FEM::FEM(const Eigen::MatrixXd &_Vf,
         const Eigen::MatrixXi &_Ff,
         const int _inner_iter,
         const Eigen::VectorXd &_young_modulus,
         const double _poisson_ratio,
         const bool _fem_symmetric):
inner_iter(_inner_iter),
young_modulus(_young_modulus),
poisson_ratio(_poisson_ratio),
mu(young_modulus.array()/2./(1.+poisson_ratio)),
lambda(young_modulus.array()*poisson_ratio/(1.+poisson_ratio)/(1.-2.*poisson_ratio)),
Vf(_Vf),
Ff(_Ff)

{
//  std::ofstream filmtr;
//  std::ofstream filmver;
  igl::per_face_normals(Vf, Ff, FNf);

  if (_fem_symmetric)
    matrix_type = -2;
  else
    matrix_type = 11;

//  filmtr.open("face_connections.txt"); (Tanay)
//  filmtr << _Ff << endl;
//  filmtr.close();

//  filmver.open("film_vertex_locations.txt", std::ofstream::app); //(Tanay)
//  filmver << _Vf << endl;
//  filmver.close();
}

void FEM::flatten_all(const Eigen::MatrixXd &Vf_uv,
                      Eigen::MatrixXd &Xflat,
                      Eigen::MatrixXd &xflat)
{
  int nf = Ff.rows();
  Xflat.resize(nf*3,2);
  xflat.resize(nf*3,2);

  for (int fi = 0; fi<nf; ++fi)
  {
    Matrix<double,3,2> q;
    q <<Vf_uv.row(Ff(fi,0)), Vf_uv.row(Ff(fi,1)), Vf_uv.row(Ff(fi,2));
    // current configuration x is the uv
    xflat.block<3, 2>(3*fi, 0) = q;

    Matrix3d P;
    P<<Vf.row(Ff(fi,0)), Vf.row(Ff(fi,1)), Vf.row(Ff(fi,2));
    Vector3d z; z<<0,0,1; //todo: is this the right direction? depends on Vf,Ff
    Vector3d n = FNf.row(fi).transpose().normalized();

    Matrix3d rot = igl::rotation_matrix_from_directions(n, z);
    P = (P.rowwise()-P.row(0))*rot.transpose();
    if( (P.col(2).cwiseAbs().maxCoeff()>1e-6) || ((P.array() != P.array()).any()) )
    {
      cerr<<"got nans in flatten_all!"<<endl;
      exit(1);
    }

    Matrix<double,3,2> p = P.topLeftCorner<3, 2>();

    // find best fitting rotation that matches rotated triangle to the 2D
    // triangle of the previous step
    Matrix2d R;
    Vector2d t;
    double scale;
    bool includeScaling = false;
    bool includeReflections = false;
    igl::procrustes(p, q, includeScaling, includeReflections, scale, R, t);

    p = (p*R).rowwise() + t.transpose();
    // undeformed equilibrium configuration X is the deformed mesh mapped to 2D
    Xflat.block<3,2>(3*fi, 0) = p;
    if((Xflat.block<3,2>(3*fi, 0).array()!=Xflat.block<3,2>(3*fi, 0).array()).any())
    {
      cerr<<"got nans in flatten_all!"<<endl;
      exit(1);
    }
  }

}
bool FEM::solve(const Eigen::VectorXi &b,
                Eigen::MatrixXd &Vf_uv, double current_step)
{
  std::ofstream Filmvert; // Get Vf (Expect to update in each step) (Tanay)
  int nf = Ff.rows(); // Ff is face triangle connectivity. Remains constant.
  int nv = Vf.rows(); // Vf is film vertex positions. It changes.

  MatrixXd Xflat, xflat;
  flatten_all(Vf_uv, Xflat, xflat);

  precompute(Xflat);

  int nc = 2*b.rows();

  // known indices in vector
  VectorXi isConstrained; isConstrained.setZero(2*nv,1);
  for (int i=0; i<b.rows(); ++i)
  {
    isConstrained[b[i]] = 1;
    isConstrained[b[i]+nv] = 1;
  }
  VectorXi ik(nc,1); ik<<b,b.array()+nv;
  VectorXi iu(2*nv-nc,1);
  int ind = 0;
  VectorXi indexInUnknown; indexInUnknown.setConstant(2*nv,1,-1);
  for (int i=0; i<2*nv; ++i)
    if (!isConstrained[i])
    {
      indexInUnknown[i] = ind;
      iu[ind++] = i;
    }




  VectorXd x;
  x.resize(2*nv,1);
  x<<Vf_uv.col(0), Vf_uv.col(1);

  double E = elastic_energy(xflat);
  // start inner iterations
  for (int iter=0; iter<inner_iter; ++iter)
  {
    if (E!=E)
    {
      cerr<<"FEM::solve: got nan in Energy"<<endl;
      exit(1);
    }
    // solve quasistatic problem. In this case, dampening
    // forces are zero (since we assume zero velocity in the settled state)
    // so we only need to use elastic forces)

    printf("Iteration %d/%d: E = %.10e\n", iter, inner_iter, E);

    // Eigen Solution

    VectorXd Dx, rhs;
    // known displacements in vector
    VectorXd Dxk; Dxk.setZero(nc,1);
    VectorXd fe;
    SparseMatrix<double> dfedx;
    elastic_forces(xflat, fe, true, dfedx, current_step);

    SparseMatrix<double> K = -dfedx;

    SparseMatrix<double> Kuu, Kuk;
    VectorXd fu;
    igl::slice(K,iu,iu,Kuu);
    igl::slice(K,iu,ik,Kuk);
    igl::slice(fe,iu,1,fu);
    rhs = fu-Kuk*Dxk;

    bool deco_success;
    if (matrix_type == 11)
    {
      Eigen::SparseLU<Eigen::SparseMatrix<double> > solverLU;
      solverLU.compute(Kuu);deco_success = solverLU.info()==Success;
      Dx = solverLU.solve(rhs);
    }
    else
    {
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solverLDLT;
      solverLDLT.compute(Kuu);deco_success = solverLDLT.info()==Success;
      Dx = solverLDLT.solve(rhs);
    }
    if(!deco_success)
    {
      cerr<<"decomposition failed!"<<endl;
      exit(1);
    }
    cerr<<"deviation --  Eigen: "<< (Kuu*Dx-rhs).cwiseAbs().maxCoeff()<<endl;



    double lambda_newton = 1;
    Eigen::VectorXd xnew = x;
    while (true)
    {
      VectorXd t = igl::slice(x, iu, 1) + lambda_newton*Dx;
      igl::slice_into(t, iu, 1, xnew);
      Vf_uv.col(0) = xnew.head(nv);
      Vf_uv.col(1) = xnew.tail(nv);
      for (int fi=0; fi<nf; ++fi)
        xflat.block<3,2>(3*fi, 0) <<Vf_uv.row(Ff(fi,0)), Vf_uv.row(Ff(fi,1)), Vf_uv.row(Ff(fi,2));

      if (elastic_energy(xflat)<E || lambda_newton <= 1e-8)
      {
        x = xnew;
        break;
      }
      lambda_newton = lambda_newton/2.;
    }


    Vf_uv.col(0) = x.head(nv);
    Vf_uv.col(1) = x.tail(nv);
    for (int fi=0; fi<nf; ++fi)
      xflat.block<3,2>(3*fi, 0) <<Vf_uv.row(Ff(fi,0)), Vf_uv.row(Ff(fi,1)), Vf_uv.row(Ff(fi,2));

    double oldE = E;
    E = elastic_energy(xflat);

    if(fabs(oldE-E)<1e-5)
      break;

  }
//  Filmvert.open("film_vertex_locations_ongoing.txt", std::ofstream::app); (Tanay)
//  Filmvert << Vf <<endl;
//  Filmvert.close();

  return true;
}

void FEM::precompute(const Eigen::MatrixXd &X)
{
  int nf = X.rows()/3;
  Bm.setZero(2*nf,2);
  W.setZero(nf,1);
  Matrix2d Dm;
  for (int fi=0; fi<nf; ++fi)
  {
    const Matrix<double, 2, 3> &X_ = X.block<3,2>(3*fi, 0).transpose();
    Dm.col(0) << X_.col(0) - X_.col(2);
    Dm.col(1) << X_.col(1) - X_.col(2);
    Bm.block<2,2>(2*fi, 0) = Dm.inverse();
    if ( (Bm.block<2,2>(2*fi, 0).array() != Bm.block<2,2>(2*fi, 0).array()). any() )
    {
      cerr<<"got nans in precompute!"<<endl;
      exit(1);
    }
    W(fi) = .5*fabs(Dm.determinant()); // triangle area
  }
}


double FEM::elastic_energy(const Eigen::MatrixXd &x)
{
  int nf = x.rows()/3;
  Matrix2d I = Matrix2d::Identity();

  double E = 0;
  for (int fi=0; fi<nf; ++fi)
  {

    Vector4d F = deformation_gradients(x, fi, false);
    Matrix2d Fmat;
    Fmat.row(0) << F.head<2>().transpose();
    Fmat.row(1) << F.tail<2>().transpose();

    Matrix2d R, S;
    igl::polar_svd(Fmat, R, S);

    Matrix2d SI = S-I;
    double PSI = mu[fi]*(SI.transpose()*SI).trace() + lambda[fi]*powl((SI).trace(),2);
    assert(PSI>=0);

    E += W(fi)*PSI;

  }
  return E;
}

Vector4d FEM::deformation_gradients(const MatrixXd &x,
                                    const int fi,
                                    bool do_derivative,
                                    Eigen::Matrix<double,4,6> &dFdx)
{
  const Matrix<double, 2, 3> &x_ = x.block<3,2>(3*fi, 0).transpose();
  const Matrix2d Bm_ = Bm.block<2,2>(2*fi, 0);

  Matrix2d Ds;
  Ds.col(0) << x_.col(0) - x_.col(2);
  Ds.col(1) << x_.col(1) - x_.col(2);
  Matrix2d Fmat = Ds*Bm_;
  //vectorised, row-major
  Eigen::Vector4d F;
  F<<Fmat.row(0).transpose(), Fmat.row(1).transpose();
  if ((F.array() != F.array()).any())
  {
    cerr<<"got nans in deformation_gradients!"<<endl;
    exit(1);
  }

  if (do_derivative)
  {
    dFdx.setZero();
    const RowVector2d &v1 = Bm_.row(0);
    const RowVector2d &v2 = Bm_.row(1);
    const RowVector2d v3 = -v1-v2;
    for (int i = 0; i<2; ++i)
    {
      dFdx(i,0)=v1(i);
      dFdx(i,1)=v2(i);
      dFdx(i,2)=v3(i);
      dFdx(2+i,3+0)=v1(i);
      dFdx(2+i,3+1)=v2(i);
      dFdx(2+i,3+2)=v3(i);
    }
  }
  return F;

}

Vector4d FEM::corotational_elasticity_stress_tensor(const MatrixXd &x,
                                                    const int fi,
                                                    bool do_derivative,
                                                    Eigen::Matrix<double,4,6> &dPdx,
                                                    double current_step)
{
 // std::ofstream strn1; //Tanay
 // std::ofstream strn2; //Tanay
  Eigen::Vector4d F;
  Matrix<double, 4, 6> dFdx;
  F = deformation_gradients(x, fi, true, dFdx);

  Matrix2d Fmat;
  Fmat.row(0) << F.head<2>().transpose();
  Fmat.row(1) << F.tail<2>().transpose();

  Matrix2d R, S;
  igl::polar_svd(Fmat, R, S);
  if (R.determinant() <0)
  {
    cerr<<"negative determinant"<<endl;
  }

  Eigen::Vector4d Rvec;
  Rvec<<R.row(0).transpose(), R.row(1).transpose();
  Matrix2d I = Matrix2d::Identity();

  Matrix2d Pmat = R*(2*mu[fi]*(S-I) + lambda[fi]*(S-I).trace()*I);
  Eigen::Vector4d P;
  P<<Pmat.row(0).transpose(), Pmat.row(1).transpose();
  if (do_derivative)
  {
    Matrix4d dRdF;
    dRdF_func(R, S, dRdF);
    Matrix<double, 4, 6> dRdx;
    dRdx = dRdF*dFdx;
    double Tr = (R.transpose()*Fmat).trace();
    Matrix<double, 1, 6> dTrdx;
    dTrdx = F.transpose()*dRdx + Rvec.transpose()*dFdx;
    dPdx = 2*mu[fi]*(dFdx-dRdx) + lambda[fi]*(Rvec*dTrdx + Tr*dRdx - I.trace()*dRdx);
  }

  std::ofstream strn1;

  if (current_step == 2)
  {
    strn1.open("strainsSI_stp2.txt", std::ofstream::app); // Repeated at each inner iteration (Tanay)
    strn1 << S-I <<endl; // strains1.txt has 2 value for each element at each inner iteration
    strn1.close();
  }
  
  if (current_step == 3)
  {
    strn1.open("strainsSI_stp3.txt", std::ofstream::app); // Repeated at each inner iteration (Tanay)
    strn1 << S-I <<endl; // strains1.txt has 2 value for each element at each inner iteration
    strn1.close();
  }

  if (current_step == 33)
  {
    strn1.open("strainsSI_stp33.txt", std::ofstream::app); // Repeated at each inner iteration (Tanay)
    strn1 << S-I <<endl; // strains1.txt has 2 value for each element at each inner iteration
    strn1.close();
  }

  if (current_step == 66)
  {
    strn1.open("strainsSI_stp66.txt", std::ofstream::app); // Repeated at each inner iteration (Tanay)
    strn1 << S-I <<endl; // strains1.txt has 2 value for each element at each inner iteration
    strn1.close();
  }

  if (current_step == 99)
  {
    strn1.open("strainsSI_stp99.txt", std::ofstream::app); // Repeated at each inner iteration (Tanay)
    strn1 << S-I <<endl; // strains1.txt has 2 value for each element at each inner iteration
    strn1.close();
  }

  this->S = S;
  this->I = I;
  return P;
}

void FEM::elastic_forces(const Eigen::MatrixXd &x,
                         Eigen::VectorXd &f,
                         bool do_derivative,
                         Eigen::SparseMatrix<double> &dfdx,
                         double current_step)
{
  int nf = Ff.rows();
  int nv = Ff.maxCoeff()+1; // the real number of vertices
  f.setZero(2*nv,1);
  Eigen::Matrix<double,4,6> dPdx;
  Vector4d P;
  vector< Triplet<double> > triplets;

  for (int fi = 0; fi<nf; ++fi)
  {
    const Matrix2d Bm_ = Bm.block<2,2>(2*fi, 0);
    P = corotational_elasticity_stress_tensor(x, fi, do_derivative, dPdx, current_step);
    Matrix2d Pmat;
    Pmat.row(0) << P.head<2>().transpose();
    Pmat.row(1) << P.tail<2>().transpose();

    Matrix2d H = -W(fi) * Pmat * Bm_.transpose();
    for (int i = 0; i<2; ++i)//indices in forces
      for (int j = 0; j<2; ++j)//indices in x,y for forces
      {
        int II = j*nv+Ff(fi,i);//row of dfdx for current corner of face
        int II1 = j*nv+Ff(fi,2);//row of dfdx for third corner of face
        f[II] += H(j,i);
        f[II1] -= H(j,i);
      }
    if (do_derivative)
    {
      Eigen::Matrix<double,4,6> dHdx;
      dHdx.topRows<2>() = -W(fi)*Bm_*dPdx.topRows<2>();
      dHdx.bottomRows<2>() = -W(fi)*Bm_*dPdx.bottomRows<2>();


      for (int i = 0; i<2; ++i)//indices in forces
      {
        for (int j = 0; j<2; ++j)//indices in x,y for forces
        {
          int ii = j*2+i;//row of dHdx
          int II = j*nv+Ff(fi,i);//row of dfdx for current corner of face
          int II1 = j*nv+Ff(fi,2);//row of dfdx for third corner of face

          for (int k = 0; k<6; ++k) // column of dHdx
          {
            int ixy = k/3;//index in x,y for x
            int ic = k%3;//index in face corner for x
            int JJ = ixy*nv+Ff(fi,ic); // column of dfdx
            triplets.push_back(Triplet<double>(II,JJ,dHdx(ii,k)));
            triplets.push_back(Triplet<double>(II1,JJ,-dHdx(ii,k)));
          }
        }
      }

    }

  }
  if (do_derivative)
  {
    dfdx.resize(2*nv,2*nv);
    dfdx.setFromTriplets(triplets.begin(), triplets.end());
  }

}
