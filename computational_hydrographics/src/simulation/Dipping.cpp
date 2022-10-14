#include "Dipping.h"

using namespace std;
using namespace Eigen;
#include <igl/arap.h>
#include <igl/barycentric_coordinates.h>
#include <igl/triangle_triangle_adjacency.h>
#include <unordered_set>
#include <igl/is_border_vertex.h>
#include <igl/per_vertex_normals.h>
#include <igl/writeOBJ.h>
#include <igl/embree/line_mesh_intersection.h>
#include <igl/barycentric_to_global.h>
#include <igl/avg_edge_length.h>
#include <igl/fit_plane.h>
#include <igl/point_in_poly.h>
#include <cmath>
#include <igl/average_onto_vertices.h>
#include <igl/viewer/Viewer.h>
#include <algorithm>
#include <igl/doublearea.h>
#include <igl/principal_curvature.h>
#include <igl/per_vertex_attribute_smoothing.h>

#include "FEM.h"
#include "utils.h"
#include <iostream>
#include <fstream>


Dipping::Dipping()
{
  current_view = 1;

  steps = 100;
  do_ARAP = false;
  resample_uv = true;
  stretch_multiplier = 0.65;
  border_fixed = 1;

  inner_iter = 10;
  young_modulus_nominal = 0.1;
  do_adaptive_young = true;
  min_young = 1e-2;
  max_young = 6;

  poisson_ratio = 0.;

  fem_symmetric = false;

  plasticity = true;

  zprojection = false;
  discoloration = true;
  discoloration_multiplier = 1;

}

Dipping::~Dipping()
{
}

void Dipping::set_ink_density_from_rgb(const Matrix<unsigned char,Dynamic,Dynamic> R,
  const Matrix<unsigned char,Dynamic,Dynamic> G,
  const Matrix<unsigned char,Dynamic,Dynamic> B)
{
  //work first with double image, then set d.IDf
  MatrixXd IDd(R.rows(), R.cols());
  for (unsigned i=0;i<IDd.rows();++i)
  {
    for (unsigned j=0;j<IDd.cols();++j)
    {
      double r = double(R(i,j))/255.0;
      double g = double(G(i,j))/255.0;
      double b = double(B(i,j))/255.0;
      double c,m,y,k;
      rgb2cmyk(r,g,b,c,m,y,k);
      IDd(i,j) = (c+m+y+k);
    }
  }
  double minId = IDd.minCoeff();
  double maxId = IDd.maxCoeff();
  IDd = (IDd.array()-minId)/(maxId-minId);

  //set the background ink density to sth constant (set to .2165 from checkerboard.png)
  double bg_id =   0.2165;
  for (unsigned i=0;i<IDd.rows();++i)
  for (unsigned j=0;j<IDd.cols();++j)
  {
    bool is_white = (R(i,j)==255 && G(i,j)==255 && B(i,j)==255);
    if (is_white)
    IDd(i,j) = bg_id;
  }

  //blur
  MatrixXd IDd_before = IDd;
  int boxHalfSize = std::max(IDd.rows(),IDd.cols())/200.;
  for (int i=0;i<IDd.rows();++i)
  for (int j=0;j<IDd.cols();++j)
  {
    int num = 0;
    double val = 0;
    for (int ii=i-boxHalfSize;ii<=i+boxHalfSize;++ii)
    for (int jj=j-boxHalfSize;jj<=j+boxHalfSize;++jj)
    {
      if (ii<0 || jj<0 || ii>=IDd.rows() || jj>=IDd.cols())
      continue;
      num++;
      val += IDd_before(ii,jj);
    }
    IDd(i,j) = val/num;
  }

  //fill in unsigned char image
  IDf.setConstant(R.rows(), R.cols(), bg_id);
  for (unsigned i=0;i<IDd.rows();++i)
  for (unsigned j=0;j<IDd.cols();++j)
  IDf(i,j) = (unsigned char) (IDd(i,j)*255.);

}
void Dipping::calculate_young_modulus_values()
{
  young_modulus_array.setConstant(Ff.rows(), 1, young_modulus_nominal);

  if (do_adaptive_young)
  {
    for (int fi=0; fi<Ff.rows(); ++fi)
    {
      //get uv at midpoint of triangle
      Eigen::Matrix<int, 3, 2> points;

      unsigned int min_x = 1e8, min_y = 1e8, max_x = 0, max_y = 0 ;
      std::vector<std::vector<unsigned int > > poly;
      for (int i =0; i<3; ++i)
      {
        std::vector<unsigned int > pt;
        pt.push_back( (unsigned int)(Vf_uv_undistorted(Ff(fi,i),0)* IDf.rows()+.5) );
        pt.push_back( (unsigned int)(Vf_uv_undistorted(Ff(fi,i),1)* IDf.cols()+.5) );
        if (min_x > pt[0])
        min_x = pt[0];
        if (min_y > pt[1])
        min_y = pt[1];
        if (max_x < pt[0])
        max_x = pt[0];
        if (max_y < pt[1])
        max_y = pt[1];
        poly.push_back(pt);
      }

      max_x = (max_x>=IDf.rows())? IDf.rows()-1: max_x;
      max_y = (max_y>=IDf.cols())? IDf.cols()-1: max_y;


      double ink_density = 0.; int num = 0;
      for (unsigned int x = min_x; x<= max_x; ++x)
      {
        for (unsigned int y = min_y; y<= max_y; ++y)
        {
          if (igl::point_in_poly(poly, x, y))
          {
            int idf = IDf(x, y);
            double val = idf;//(idf<0)?256.-fabs(idf):idf;
            val = val/255.;
            ink_density += val;
            num++;
          }
        }
      }
      ink_density/= num;

      young_modulus_array[fi] = min_young + ink_density*(max_young-min_young);
    }
  }
}

void Dipping::init_texture_mesh(double squaresize, int u_count, int v_count) //made (10.75,20,20) (Tanay)
{
  int cols = v_count*15; //default: 15.0 (Tanay)
  int rows = u_count*15; //default: 15.0 (Tanay)

  auto ig = [&] (int i, int j) {return j*rows + i;};

  Vf.resize(rows*cols,3);
  for (unsigned i=0;i<rows;++i)
  for (unsigned j=0;j<cols;++j)
  Vf.row(ig(i,j)) << double(i)/(rows-1),double(j)/(cols-1),0;

  Vf_uv = Vf.block(0,0,Vf.rows(),2);
  Vf_uv_undistorted = Vf.block(0,0,Vf.rows(),2);

  scale_u = squaresize * u_count;
  scale_v = squaresize * v_count;

  // Scale to the correct size
  Vf.col(0) = Vf.col(0).array() * scale_u;
  Vf.col(1) = Vf.col(1).array() * scale_v;

  // Center the grid
  for (unsigned i=0; i<2;++i)
  Vf.col(i) = Vf.col(i).array() - (Vf.col(i).maxCoeff()+Vf.col(i).minCoeff())/2.0;

  Ff.resize((cols-1)*(rows-1)*2,3);
  int counter = 0;
  for (unsigned i=0;i<rows-1;++i)
  {
    for (unsigned j=0;j<cols-1;++j)
    {
      Ff.row(counter++) << ig(i,j), ig(i+1,j), ig(i+1,j+1);
      Ff.row(counter++) << ig(i+1,j+1), ig(i,j+1), ig(i,j);
    }
  }

  // Compute and store the original UV coordinate and triangle (in UV space) of all the
  // vertices of the texture mesh
  Vf_uv_ori = Vf_uv;
  Vf_uv_fid.resize(Vf_uv.rows());
  for(unsigned i=0; i<Ff.rows(); ++i)
  {
    if (Ff(i,0) >=0) Vf_uv_fid(Ff(i,0)) = i;
    if (Ff(i,1) >=0) Vf_uv_fid(Ff(i,1)) = i;
    if (Ff(i,2) >=0) Vf_uv_fid(Ff(i,2)) = i;
  }

  igl::triangle_triangle_adjacency(Vf,Ff,TTf,TTif);

  border = igl::is_border_vertex(Vf,Ff);

  Vf_ori = Vf;

  calculate_young_modulus_values();
}

bool in_triangle_3d(const Vector3d& v1, const Vector3d& v2, const Vector3d& v3, const Vector3d& p, Vector3d& bc)
{
  Vector3d a3 = (v2-v1).cross(v3-v1);
  a3.normalize();

  Vector3d a1 = (v2-v1).normalized();
  Vector3d a2 = a3.cross(a1).normalized();

  // Project

  Vector2d p_2(p.dot(a1), p.dot(a2));
  Vector2d v1_2(v1.dot(a1), v1.dot(a2));
  Vector2d v2_2(v2.dot(a1), v2.dot(a2));
  Vector2d v3_2(v3.dot(a1), v3.dot(a2));

  Vector3d z1,z2,z3,z4;
  z1 << v1_2(0), v1_2(1), 1;
  z2 << v2_2(0), v2_2(1), 1;
  z3 << v3_2(0), v3_2(1), 1;
  MatrixXd A(3,3);
  A << z1, z2, z3;
  Vector3d b;
  b << p_2(0), p_2(1),1;

  Vector3d t = A.colPivHouseholderQr().solve(b);

  bc = t;
  return (t.minCoeff() > -0.1);

}

void Dipping::init()
{
  // Project reference mesh for comparisons
  if (Vq.rows() > 0)
  {
    if (false)
    {
      MatrixXd N = MatrixXd::Zero(Vq.rows(),3);
      for (unsigned i=0;i<Fq.rows();++i)
      {
        MatrixXd V_quad(4,3);
        V_quad << Vq.row(Fq(i,0)),Vq.row(Fq(i,1)),Vq.row(Fq(i,2)),Vq.row(Fq(i,3));
        RowVector3d n;
        RowVector3d c;

        igl::fit_plane(V_quad,n,c);
        N.row(i) = n;
      }
      MatrixXd Vq_temp = Vq + (N.array()*0.001).matrix();

      Vq_bc = igl::embree::line_mesh_intersection(Vq_temp,N,V,F);

    }
    else
    {
      Vq_bc = MatrixXd::Zero(Vq.rows(),4);

      // Project the point on each triangle and check bc
      for(unsigned z=0; z<Vq.rows();++z)
      {
        double closest = 100000000;

        Vector3d bc_temp;
        for(unsigned i=0; i<F.rows();++i)
        {
          if (in_triangle_3d(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),Vq.row(z),bc_temp))
          {
            RowVectorXd p_temp = V.row(F(i,0))*bc_temp(0) + V.row(F(i,1))*bc_temp(1) + V.row(F(i,2))*bc_temp(2);
            if ((p_temp - Vq.row(z)).norm() < closest)
            {
              closest = (p_temp - Vq.row(z)).norm();
              Vq_bc.row(z) << i, bc_temp.transpose();
            }
          }
        }
      }
    }
  }
}

void Dipping::init_simulation()
{
  init_texture_mesh(10.75, 20, 20); // 10 (originally 25.25, 8, 11) (Tanay)

  // unfreeze all vertices
  fixed = VectorXi::Zero(Vf.rows());
  new_fixed = VectorXi::Zero(Vf.rows());

  // Move the object to touch the water
  if (Vq.rows() > 0)
  Vq.col(2) = Vq.col(2).array() - V.col(2).minCoeff();

  V.col(2) = V.col(2).array() - V.col(2).minCoeff();

  // Dipping speed requires 100 iterations
  dt = (V.col(2).maxCoeff() - V.col(2).minCoeff())/steps;

  // Init uv
  //Olga: these seem to be repeated here. to remove?
  Vf_uv = Vf.block(0,0,Vf.rows(),2);
  Vf_uv_ori = Vf_uv;

  // Reset the steps
  current_step = 0;

  // Copy the colors for Rf
  Rf = Rf_ori;
  Gf = Gf_ori;
  Bf = Bf_ori;

  // Clean the final texture
  R_final.setZero(R.rows(),R.cols());
  G_final.setZero(R.rows(),R.cols());
  B_final.setZero(R.rows(),R.cols());

  std::ofstream filmtrinit;
  std::ofstream filmverinit;

  filmtrinit.open("face_connections.txt");
  filmtrinit << Ff << endl;
  filmtrinit.close();

  filmverinit.open("film_vertex_locations_stp0.txt", std::ofstream::app);
  filmverinit << Vf << endl;
  filmverinit.close();

}

bool Dipping::simulate()
{
  current_step++;

  cerr<<"--- Current step: "<<current_step<<"--- "<<endl;
  if (current_step > steps)
  return false;

  // Move the object down
  // Every vertex of the grid contained in triangles with z<=0 are fixed
  // and move down
  move_object_and_film_down_project();

  if (!zprojection)
  {
    // Solve ARAP on the deformed film with fixed UVs
    if (!deform())
    return false;
  }

  // Resample UV mesh
  if (resample_uv)
  resample();

  return true;
}

void Dipping::move_object_and_film_down_project()
{
  // Move down the object
  V.col(2) = V.col(2).array() - dt;

  // Project Vf onto V
  MatrixXd N(Vf.rows(),3);
  for (unsigned i=0;i<Vf.rows();++i)
  N.row(i) << 0,0,1;

  MatrixXd BC = igl::embree::line_mesh_intersection(Vf,N,V,F);
  MatrixXd Vftemp = igl::barycentric_to_global(V,F,BC);

  for (unsigned i=0; i<Vftemp.rows();++i)
  {
    if (fixed(i))
    Vf(i,2) = Vf(i,2) - dt;
    else
    {
      if (Vftemp(i,2) < 0)
      {
        Vf(i,2) = max(Vftemp(i,2),Vf(i,2)-dt);
        new_fixed(i) = 1;
      }
    }
  }
}

bool Dipping::deform()
{
  MatrixXd debuguv = Vf_uv;

  bool success;
  if (do_ARAP)
  success = deformARAP();
  else
  success = deformFEM();

  for (unsigned i=0;i<fixed.size();++i)
  {
    if (fixed(i) != 0)
    {
      assert((Vf_uv.row(i)-debuguv.row(i)).norm() < 1e-10);
    }
  }

  for (unsigned i=0;i<new_fixed.size();++i)
  if (new_fixed(i) != 0)
  fixed(i) = 1;

  new_fixed = VectorXi::Zero(Vf.rows());

  return success;

}
bool Dipping::deformFEM()
{
  int nconst = fixed.sum();
  if (border_fixed)
  for (unsigned i=0; i<border.size(); ++i)
  if (border[i])
  nconst += 1;
  Eigen::VectorXi b  = Eigen::VectorXi::Zero(nconst);

  int count = 0;
  for (unsigned i=0; i<fixed.size(); ++i)
  if (fixed(i) != 0 || (border[i] && border_fixed))
  {
    b(count) = i;count++;
  }


  MatrixXd Vf_temp = Vf;
  Vf_temp.col(2) = Vf_temp.col(2).array() * stretch_multiplier;
  
  FEM fem(Vf_temp, Ff, inner_iter, young_modulus_array, poisson_ratio, fem_symmetric);

  bool retval = fem.solve(b, Vf_uv, current_step);

  std::ofstream filmver;

  if (current_step == 2)
  {
    filmver.open("film_vertex_locations_stp2.txt", std::ofstream::app);
    filmver << Vf << endl;
    filmver.close();
  }

  if (current_step == 3)
  {
  	filmver.open("film_vertex_locations_stp3.txt", std::ofstream::app);
    filmver << Vf << endl;
    filmver.close();
  }

  if (current_step == 33)
  {
    filmver.open("film_vertex_locations_stp33.txt", std::ofstream::app);
    filmver << Vf << endl;
    filmver.close();
  }

  if (current_step == 66)
  {
    filmver.open("film_vertex_locations_stp66.txt", std::ofstream::app);
    filmver << Vf << endl;
    filmver.close();
  }

  if (current_step == 99)
  {
    filmver.open("film_vertex_locations_stp99.txt", std::ofstream::app);
    filmver << Vf << endl;
    filmver.close();
  }

  return retval;
}

bool Dipping::deformARAP()
{
  // Solve ARAP with fixed vertices (ARAP: As rigid as possible) (Tanay)

  // Add dynamic regularization to avoid to specify boundary conditions
  igl::ARAPData arap_data;
  arap_data.with_dynamics = false;

  // count the fixed vertices and add them as constraints
  int nconst = fixed.sum();

  if (border_fixed)
  for (unsigned i=0; i<border.size(); ++i)
  if (border[i])
  nconst += 1;

  Eigen::VectorXi b  = Eigen::VectorXi::Zero(nconst);
  Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(nconst,2);

  int count = 0;
  for (unsigned i=0; i<fixed.size(); ++i)
  {
    if (fixed(i) != 0 || (border[i] && border_fixed))
    {
      b(count) = i;
      bc.row(count) = Vf_uv.block(i,0,1,2);
      ++count;
    }
  }

  // Initialize ARAP
  arap_data.max_iter = 10;

  // Solve ARAP

  MatrixXd Vf_temp;

  if (plasticity)
  {
    Vf_temp = Vf;
    Vf_temp.col(2) = Vf_temp.col(2).array() * stretch_multiplier;
  }
  else
  {
    Vf_temp = Vf_ori;
  }

  arap_precomputation(Vf_temp,Ff,2,b,arap_data);
  arap_solve(bc,arap_data,Vf_uv);
  return true;
}

void resample_uv_aux(const MatrixXi& TT, const MatrixXd& UV, const MatrixXd& V, const MatrixXi& F, const RowVector2d& p, int face_hint, Vector3d& res, int& new_hint)
{
  if (false)
  {
    for (unsigned i=0; i<F.rows();++i)
    {
      Vector3d v1,v2,v3,v4;
      v1 << UV(F(i,0),0), UV(F(i,0),1), 1;
      v2 << UV(F(i,1),0), UV(F(i,1),1), 1;
      v3 << UV(F(i,2),0), UV(F(i,2),1), 1;
      MatrixXd A(3,3);
      A << v1, v2, v3;
      Vector3d b;
      b << p(0), p(1),1;

      Vector3d t = A.colPivHouseholderQr().solve(b);

      if (t.minCoeff() > -0.0001)
      {
        new_hint = i;
        res = V.row(F(i,0)) * t(0) + V.row(F(i,1)) * t(1) + V.row(F(i,2)) * t(2);

        return;
      }
    }
    assert(false);
  }
  else
  {
    double previous_min = -100;
    int previous_hint = -1;
    Vector3d previous_res;

    std::vector<int> Q;
    Q.reserve(25);

    std::unordered_set<int> visited;
    size_t iter = 0;

    Q.push_back(face_hint);
    visited.insert(face_hint);
    while(iter<Q.size())
    {
      int fid = Q[iter++];

      Vector3d v1,v2,v3,v4;
      v1 << UV(F(fid,0),0), UV(F(fid,0),1), 1;
      v2 << UV(F(fid,1),0), UV(F(fid,1),1), 1;
      v3 << UV(F(fid,2),0), UV(F(fid,2),1), 1;
      MatrixXd A(3,3);
      A << v1, v2, v3;
      Vector3d b;
      b << p(0), p(1),1;

      Vector3d t = A.colPivHouseholderQr().solve(b);

      if (t.minCoeff() > -0.001)
      {
        new_hint = fid;
        res = V.row(F(fid,0)) * t(0) + V.row(F(fid,1)) * t(1) + V.row(F(fid,2)) * t(2);
        return;
      }
      else if (t.minCoeff() > previous_min)
      {
        previous_min = t.minCoeff();
        previous_hint = fid;
        previous_res = V.row(F(fid,0)) * t(0) + V.row(F(fid,1)) * t(1) + V.row(F(fid,2)) * t(2);
      }

      for (int i = 0; i<3; ++i)
      {
        int next = TT(fid,i);
        if(next >=0 &&  visited.insert(next).second)
        Q.push_back(next);
      }
    }

    if (previous_hint != -1)
    {
      new_hint = previous_hint;
      res = previous_res;
      return;
    }

    cerr << "Resampling failed, debug me!" << endl;
    assert(false);
  }
}

void Dipping::resample()
{
  MatrixXd Vf_new = Vf;
  int limit = Vf.rows();
  for (unsigned i=0;i<limit;++i)
  {
    Vector3d t;
    RowVector2d p;
    p << Vf_uv_ori(i,0), Vf_uv_ori(i,1);
    int new_hint;
    resample_uv_aux(TTf, Vf_uv, Vf, Ff, p, Vf_uv_fid(i),t,new_hint);

    Vf_uv_fid(i) = new_hint;
    Vf_new.row(i) = t;
  }
  // Replace the old mesh
  Vf = Vf_new;
  Vf_uv = Vf_uv_ori;

}

void Dipping::invert_deformation()
{
  // TODO this does not work without resample UV

  // project the original mesh onto the deformed film
  MatrixXd N_;
  igl::per_vertex_normals(V, F, N_);

  MatrixXd N;
  igl::per_vertex_attribute_smoothing(N_, F, N);

  for (int vi = 0; vi<V.rows(); ++vi)
  {
    if (N.row(vi).sum() != N.row(vi).sum())
    {
      cerr<<"invert_deformation: nans in normals"<<endl;
      N.row(vi)<< 0,0,1;
    }
  }
  MatrixXd BC = igl::embree::line_mesh_intersection(V,N,Vf,Ff);

  MatrixXd V_hit = igl::barycentric_to_global(Vf,Ff,BC);
  VectorXd d = (V_hit - V).rowwise().norm();
  VectorXi mark = (d.array()>2.).cast<int>();

  // Convert to global coordinates
  if (false)
  {
    V_uv_flat = igl::barycentric_to_global(Vf,Ff,BC);
    F_uv_flat = F;
    F_uv_flat_valid.setOnes(F_uv_flat.rows(),1);
  }
  else
  {

    MatrixXd Vftemp = MatrixXd::Zero(Vf_uv_undistorted.rows(), 3);
    if (resample_uv)
    {
      Vftemp.col(0) = Vf_uv_undistorted.col(0);
      Vftemp.col(1) = Vf_uv_undistorted.col(1);
    }
    else
    {
      MatrixXd tempuv = Vf_uv;
      for (unsigned i=0; i<2;++i)
      tempuv.col(i) = tempuv.col(i).array() - tempuv.col(i).minCoeff();

      // Scale to the correct size
      tempuv.col(0) = tempuv.col(0).array() / scale_u;
      tempuv.col(1) = tempuv.col(1).array() / scale_v;

      Vftemp.col(0) = tempuv.col(0);
      Vftemp.col(1) = tempuv.col(1);
    }

    MatrixXd temp = igl::barycentric_to_global(Vftemp,Ff,BC);

    // Center the grid with respect to the plane center
    for (unsigned i=0; i<2;++i)
    temp.col(i) = temp.col(i).array() - (Vf_uv_undistorted.col(i).maxCoeff()+Vf_uv_undistorted.col(i).minCoeff())/2.0;

    // Scale to the correct size
    temp.col(0) = temp.col(0).array() * scale_u;
    temp.col(1) = temp.col(1).array() * scale_v;

    V_uv_flat = MatrixXd::Zero(temp.rows(),3);
    V_uv_flat.col(0) = temp.col(0);
    V_uv_flat.col(1) = temp.col(1);

    F_uv_flat = F;
    F_uv_flat_valid.setOnes(F_uv_flat.rows(),1);

  }

  // Clean up degenerate faces due to projections errors
  double avg = igl::avg_edge_length(V,F);

  for (unsigned i=0; i<F_uv_flat.rows();++i)
  {
    double e1 = (V_uv_flat.row(F_uv_flat(i,0)) - V_uv_flat.row(F_uv_flat(i,1))).norm();
    double e2 = (V_uv_flat.row(F_uv_flat(i,1)) - V_uv_flat.row(F_uv_flat(i,2))).norm();
    double e3 = (V_uv_flat.row(F_uv_flat(i,2)) - V_uv_flat.row(F_uv_flat(i,0))).norm();
    if (e1 > 10*avg || e2 > 10*avg || e3 > 10*avg)
    {
      F_uv_flat.row(i) << 0,0,0;
      F_uv_flat_valid[i] = 0;
    }
  }

  // Filter everything that projected too far
  for(unsigned i=0; i<F_uv_flat.rows(); ++i)
  {
    if ( mark(F_uv_flat(i,0)) >0 || mark(F_uv_flat(i,1))>0 || mark(F_uv_flat(i,2))>0 )
    {
      F_uv_flat.row(i) << 0,0,0;
      F_uv_flat_valid[i] = 0;
    }
  }

  // Filter everything that is flipped
  for(unsigned i=0; i<F_uv_flat.rows(); ++i)
  {
    Vector3d p0 = V_uv_flat.row(F_uv_flat(i,0));
    Vector3d p1 = V_uv_flat.row(F_uv_flat(i,1));
    Vector3d p2 = V_uv_flat.row(F_uv_flat(i,2));

    Vector3d t = (p1-p0).cross(p2-p0);

    if (t(2) > 0)
    {
      F_uv_flat.row(i) << 0,0,0;
      F_uv_flat_valid[i] = 0;
    }
  }

  // If available, deform also the reference quad mesh
  if (Vq.rows() > 0)
  {
    Vq_def = MatrixXd::Zero(Vq.rows(), 3);

    for (unsigned i=0;i<Vq.rows();++i)
    Vq_def.row(i) = V_uv_flat.row(F_uv_flat(Vq_bc(i,0),0))*Vq_bc(i,1) + V_uv_flat.row(F_uv_flat(Vq_bc(i,0),1))*Vq_bc(i,2) + V_uv_flat.row(F_uv_flat(Vq_bc(i,0),2))*Vq_bc(i,3);
  }

}

double Dipping::measure_error()
{
  double ideal_edge = 25.25/2.0;
  double ideal_diag = sqrt(ideal_edge*ideal_edge + ideal_edge*ideal_edge);

  double error = 0;

  for (unsigned i=0;i<Fq.rows();++i)
  {
    // for every quad, emasure the edges
    error += 0.5*pow(((Vq_def.row(Fq(i,0)) - Vq_def.row(Fq(i,1))).norm() - ideal_edge),2);
    error += 0.5*pow(((Vq_def.row(Fq(i,1)) - Vq_def.row(Fq(i,2))).norm() - ideal_edge),2);
    error += 0.5*pow(((Vq_def.row(Fq(i,2)) - Vq_def.row(Fq(i,3))).norm() - ideal_edge),2);
    error += 0.5*pow(((Vq_def.row(Fq(i,3)) - Vq_def.row(Fq(i,0))).norm() - ideal_edge),2);

    error += pow(((Vq_def.row(Fq(i,0)) - Vq_def.row(Fq(i,2))).norm() - ideal_diag),2);
    error += pow(((Vq_def.row(Fq(i,1)) - Vq_def.row(Fq(i,3))).norm() - ideal_diag),2);
  }

  return sqrt(error);
}


void Dipping::update_discoloration_forward()
{
  if (Vf_ori.rows() == 0)
  return;


  if (!discoloration)
  {
    Rf = Rf_ori;
    Gf = Gf_ori;
    Bf = Bf_ori;
    return;
  }

  // Prepare the render buffer
  Rd.resize(Rf_ori.rows(),Rf_ori.cols());
  Gd.resize(Rf_ori.rows(),Rf_ori.cols());
  Bd.resize(Rf_ori.rows(),Rf_ori.cols());

  // Compute per triangle area ratios
  MatrixXd A_ori,A, D, D_vert;

  igl::doublearea(Vf_ori, Ff, A_ori);
  igl::doublearea(Vf, Ff, A);

  D = A_ori.array()/A.array();

  // Interpolate on vertices
  igl::average_onto_vertices(Vf_ori,Ff, D, D_vert);

  // Normalize between 0 and 1
  for(unsigned i=0;i<D_vert.size();++i)
  {
    assert(D_vert(i) > 0);
    D_vert(i) = min(D_vert(i),1.);
  }

  // Render into the texture buffer
  igl::viewer::OpenGL_state opengl;
  igl::viewer::ViewerCore core;
  igl::viewer::ViewerData data;

  // Init opengl shaders and buffers
  opengl.init();

  // Init ViewerCore
  core.init();

  // Set rendering options
  core.orthographic = true;
  core.show_lines = false;
  core.show_texture = false;
  core.lighting_factor = 0;
  core.depth_test = false;
  core.invert_normals = true;

  MatrixXd Vf_scaled = Vf_ori;

  for(unsigned i=0;i<3;++i)
  Vf_scaled.col(i) = Vf_scaled.col(i).array() - Vf_ori.col(i).minCoeff();
  for(unsigned i=0;i<2;++i)
  Vf_scaled.col(i) = Vf_scaled.col(i).array() / (Vf_ori.col(i).maxCoeff() - Vf_ori.col(i).minCoeff());
  for(unsigned i=0;i<2;++i)
  Vf_scaled.col(i) = (Vf_scaled.col(i).array() *2) - 1;

  Vf_scaled.col(2) = Vf_scaled.col(2).array() * 0;

  data.set_mesh(Vf_scaled, Ff);

  MatrixXd C(Vf_ori.rows(),3);
  for (unsigned i=0;i<C.rows();++i)
  C.row(i) << D_vert(i),D_vert(i),D_vert(i);

  data.set_colors(C);

  C_vert = C;

  core.proj = Matrix4f::Identity();
  core.model = Matrix4f::Identity();
  core.view = Matrix4f::Identity();

  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Ad(Rd.rows(),Rd.cols());

  core.background_color << 0,0,0;

  // Do the rendering in a buffer
  core.draw_buffer(data,opengl,false,Rd,Gd,Bd,Ad);

  //  texture_to_png("tmp3.png", Rd, Gd, Bd);

  //  texture_to_png("tmp6.png", Rf_ori, Gf_ori, Bf_ori);

  //Prepare the output
  Rf.resize(Rf_ori.rows(),Rf_ori.cols());
  Gf.resize(Rf_ori.rows(),Rf_ori.cols());
  Bf.resize(Rf_ori.rows(),Rf_ori.cols());

  // Blend with the original colors
  for (unsigned i=0;i<Rd.rows();++i)
  {
    for (unsigned j=0;j<Rd.cols();++j)
    {
      double w = double(Rd(i,j))/255.0;

      double r = double(Rf_ori(i,j))/255.0;
      double g = double(Gf_ori(i,j))/255.0;
      double b = double(Bf_ori(i,j))/255.0;
      double c,m,y,k;
      rgb2cmyk(r,g,b,c,m,y,k);
      c = w*c;
      m = w*m;
      y = w*y;
      k = w*k;
      cmyk2rgb(c,m,y,k,r,g,b);
      Rf(i,j) = r*255;
      Gf(i,j) = g*255;
      Bf(i,j) = b*255;
    }
  }
  //  texture_to_png("tmp4.png", Rf, Gf, Bf);

}

void Dipping::update_discoloration_inverse_to_film()
{
  if (Vf_ori.rows() == 0)
  return;


  if (!discoloration)
  {
    return;
  }


  // Prepare the render buffer
  Rd.resize(Rf_ori.rows(),Rf_ori.cols());
  Gd.resize(Rf_ori.rows(),Rf_ori.cols());
  Bd.resize(Rf_ori.rows(),Rf_ori.cols());

  // Compute per triangle area ratios
  MatrixXd A_ori,A, D, D_vert;

  igl::doublearea(Vf_ori, Ff, A_ori);
  igl::doublearea(Vf, Ff, A);

  D = A_ori.array()/A.array();

  // Interpolate on vertices
  igl::average_onto_vertices(Vf_ori,Ff, D, D_vert);

  // Normalize between 0 and 1
  for(unsigned i=0;i<D_vert.size();++i)
  {
    assert(D_vert(i) > 0);
    D_vert(i) = min(D_vert(i),1.);
  }

  // Render into the texture buffer
  igl::viewer::OpenGL_state opengl;
  igl::viewer::ViewerCore core;
  igl::viewer::ViewerData data;

  // Init opengl shaders and buffers
  opengl.init();

  // Init ViewerCore
  core.init();

  // Set rendering options
  core.orthographic = true;
  core.show_lines = false;
  core.show_texture = false;
  core.lighting_factor = 0;
  core.depth_test = false;
  core.invert_normals = true;

  MatrixXd Vf_scaled = Vf_ori;

  for(unsigned i=0;i<3;++i)
  Vf_scaled.col(i) = Vf_scaled.col(i).array() - Vf_ori.col(i).minCoeff();
  for(unsigned i=0;i<2;++i)
  Vf_scaled.col(i) = Vf_scaled.col(i).array() / (Vf_ori.col(i).maxCoeff() - Vf_ori.col(i).minCoeff());
  for(unsigned i=0;i<2;++i)
  Vf_scaled.col(i) = (Vf_scaled.col(i).array() *2) - 1;

  Vf_scaled.col(2) = Vf_scaled.col(2).array() * 0;

  data.set_mesh(Vf_scaled, Ff);

  MatrixXd C(Vf_ori.rows(),3);
  for (unsigned i=0;i<C.rows();++i)
  C.row(i) << D_vert(i),D_vert(i),D_vert(i);

  data.set_colors(C);

  C_vert = C;

  core.proj = Matrix4f::Identity();
  core.model = Matrix4f::Identity();
  core.view = Matrix4f::Identity();

  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Ad(Rd.rows(),Rd.cols());

  core.background_color << 0,0,0;

  // Do the rendering in a buffer
  core.draw_buffer(data,opengl,false,Rd,Gd,Bd,Ad);

  //  texture_to_png("tmp3.png", Rd, Gd, Bd);

  //  texture_to_png("tmp6.png", Rf_ori, Gf_ori, Bf_ori);

  //Prepare the output
  Rf.resize(Rf_ori.rows(),Rf_ori.cols());
  Gf.resize(Rf_ori.rows(),Rf_ori.cols());
  Bf.resize(Rf_ori.rows(),Rf_ori.cols());

  // Blend with the original colors
  for (unsigned i=0;i<Rd.rows();++i)
  {
    for (unsigned j=0;j<Rd.cols();++j)
    {
      double w = double(Rd(i,j))/255.0;
      if (fabs(w)<1e-6)
      w = 1.0/255.0;

      double r = double(Rf_ori(i,j))/255.0;
      double g = double(Gf_ori(i,j))/255.0;
      double b = double(Bf_ori(i,j))/255.0;
      double c,m,y,k;
      rgb2cmyk(r,g,b,c,m,y,k);

      c = std::max<double>(0,std::min<double>(1.,c/w));
      m = std::max<double>(0,std::min<double>(1.,m/w));
      y = std::max<double>(0,std::min<double>(1.,y/w));
      k = std::max<double>(0,std::min<double>(1.,k/w));
      cmyk2rgb(c,m,y,k,r,g,b);

      Rf(i,j) = r*255;
      Gf(i,j) = g*255;
      Bf(i,j) = b*255;
    }
  }

}

void Dipping::update_discoloration_inverse()
{
  if (V_uv_flat.rows() == 0)
  return;

  if (!discoloration)
  {
    R_final = R;
    G_final = G;
    B_final = B;
    return;
  }
  // Prepare buffer
  //Olga: these zeros don't really matter cause the background color is set from the viewer below, but always good to initialize:)
  Rd.setZero(R.rows(),R.cols());
  Gd.setZero(R.rows(),R.cols());
  Bd.setZero(R.rows(),R.cols());

  // Compute per triangle area ratios
  MatrixXd A_ori,A, D, D_vert;

  igl::doublearea(V, F, A_ori);

  //Olga: the invalid faces should get zero distortion, so use F_uv_flat where those faces have zero area
  igl::doublearea(V_uv_flat, F_uv_flat, A);

  D = A.array()/A_ori.array();

  // Interpolate on vertices
  //  igl::average_onto_vertices(V_uv,F_uv, D, D_vert);
  //Olga: only interpolate from valid faces, otherwise we get junk (or zeros depending on what is done in doublearea above)
  D_vert.setZero(V_uv.rows(),1);
  VectorXi COUNT; COUNT.setZero(V_uv.rows(),1);
  for (int i = 0; i <F_uv.rows(); ++i)
  {
    if (!F_uv_flat_valid[i])
    continue;
    for (int j = 0; j<3; ++j)
    {
      D_vert.row(F_uv(i,j)) += D.row(i);
      COUNT[F_uv(i,j)] ++;
    }
  }

  double min_D = 1e8;
  double max_D = 0;
  for (int i = 0; i <V_uv.rows(); ++i)
  {
    if (COUNT[i]>0)
    {
      D_vert.row(i) /= COUNT[i];
      if (min_D>D_vert(i,0))
      min_D = D_vert(i,0);
      if (max_D<D_vert(i,0))
      max_D = D_vert(i,0);
    }
  }

  // do not normalize, instead cap between 0 and 1. This also helps to cap any effects from remaining invalid faces
  //  incorrectly rendered before.
  for(unsigned i=0;i<D_vert.size();++i)
  {
    assert(D_vert(i) >= 0);
    D_vert(i) = min(D_vert(i),1.);
  }
  // for option 2, need to set those cause they are used below
  min_D = 0.; max_D = 1.;

  // Render into the texture buffer
  igl::viewer::OpenGL_state opengl;
  igl::viewer::ViewerCore core;
  igl::viewer::ViewerData data;

  // Init opengl shaders and buffers
  opengl.init();

  // Init ViewerCore
  core.init();

  // Set rendering options
  core.orthographic = true;
  core.show_lines = false;
  core.show_texture = false;
  core.lighting_factor = 0;
  core.depth_test = false;
  core.invert_normals = true;
  //Olga: background of the image will be white
  core.background_color<<1.,1.,1.;

  MatrixXd V_temp = MatrixXd::Zero(V_uv.rows(),V_uv.cols());

  for(unsigned i=0;i<2;++i)
  V_temp.col(i) = (V_uv.col(i).array() *2) - 1;

  //  //Olga
  //  data.set_mesh(V_temp, F_uv);
  // Olga: Do not render invalid faces, otherwise we get z-fighting and artifacts
  int nvalid = F_uv_flat_valid.sum();
  MatrixXi valid_F_uv(nvalid,3);
  int ind = 0;
  for (int i =0;i<F_uv.rows(); ++i)
  if(F_uv_flat_valid[i])
  valid_F_uv.row(ind++) = F_uv.row(i);
  data.set_mesh(V_temp, valid_F_uv);

  MatrixXd C(V_uv.rows(),3);
  for (unsigned i=0;i<C.rows();++i)
  //Olga: no need to render a funky image here, we're only using the red channel below
  //    C.row(i) << D_vert(i),0,1-D_vert(i);
  C.row(i) << D_vert(i),D_vert(i),D_vert(i);


  data.set_colors(C);

  //Olga: remove this, not used
  C_vert = C;

  core.proj = Matrix4f::Identity();
  core.model = Matrix4f::Identity();
  core.view = Matrix4f::Identity();

  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Ad(Rd.rows(),Rd.cols());
  // Do the rendering in a buffer
  core.draw_buffer(data,opengl,false,Rd,Gd,Bd,Ad);

  //texture_to_png("tmp1.png", Rd, Gd, Bd);

  // Enhance the original colors to compensate
  for (unsigned i=0;i<Rd.rows();++i)
  {
    for (unsigned j=0;j<Rd.cols();++j)
    {
      double w;

      double val = double(Rd(i,j))/255.0;
      w = min_D + (max_D-min_D)*val;
      if (fabs(w)<1e-6)
      w = 1.0/255.0;

      double r = double(R(i,j))/255.0;
      double g = double(G(i,j))/255.0;
      double b = double(B(i,j))/255.0;
      double c,m,y,k;
      rgb2cmyk(r,g,b,c,m,y,k);

      c = c/w;
      m = m/w;
      y = y/w;
      k = k/w;
      c = max(min(c,1.),0.);
      m = max(min(m,1.),0.);
      y = max(min(y,1.),0.);
      k = max(min(k,1.),0.);

      cmyk2rgb(c,m,y,k,r,g,b);

      R_final(i,j) = r*255;
      G_final(i,j) = g*255;
      B_final(i,j) = b*255;

    }
  }

  //texture_to_png("tmp2.png", R_final, G_final, B_final);


}


void Dipping::fill_black(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B)
  {
    bool done = false;

    while (!done)
    {
      done = true;

      for(int i=0;i<R.rows();++i)
      {
        for(int j=0;j<R.cols();++j)
        {
          if (R(i,j) != 0 && G(i,j) != 0 && B(i,j) != 0)
          {
            for(int ii=-1;ii<=1;++ii)
            {
              for(int jj=-1;jj<=1;++jj)
              {
                int ti = max(min(i+ii,int(R.rows())-1),0);
                int tj = max(min(j+jj,int(R.cols())-1),0);
                if (R(ti,tj) == 0 && G(ti,tj) == 0 && B(ti,tj) == 0)
                {
                  R(i+ii,j+jj) = R(i,j);
                  G(i+ii,j+jj) = G(i,j);
                  B(i+ii,j+jj) = B(i,j);
                  done = false;
                }
              }
            }
          }
        }
      }
    }
  }
