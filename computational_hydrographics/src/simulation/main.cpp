#if 1
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include "Dipping.h"
#include "Dipping_serialize.h"

using namespace std;
using namespace Eigen;

Dipping d;
igl::viewer::Viewer viewer;
std::string line_search_dir("line_search_results_A/");

//#define TIME
#ifdef TIME
#include <igl/Timer.h>
igl::Timer timer;
#endif
#include <igl/matlab_format.h>
#include <igl/jet.h>
#include <igl/list_to_matrix.h>
#include <igl/file_dialog_save.h>
#include <igl/file_dialog_open.h>
#include <igl/png/render_to_png.h>
#include <iostream>
#include "utils.h"
#include <fstream>


void save_image(const std::string &img_filename)
{
  igl::viewer::OpenGL_state opengl;
  igl::viewer::ViewerCore core;
  igl::viewer::ViewerData data;

  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A;

  // Init opengl shaders and buffers
  opengl.init();

  // Init ViewerCore
  core.init();

  // Set rendering options
  core.orthographic = true;
  core.show_lines = false;
  core.show_texture = true;
  core.lighting_factor = 0;
  core.depth_test = false;
  core.invert_normals = true;

  // Prepare the local buffer
  RowVector3d max_size = d.Vf_ori.colwise().maxCoeff()-d.Vf_ori.colwise().minCoeff();
  double mult = 300./25.4;
  // Prepare the local buffer
  R.resize(floor(max_size(0)*mult), floor(max_size(1)*mult));
  G.resize(floor(max_size(0)*mult), floor(max_size(1)*mult));
  B.resize(floor(max_size(0)*mult), floor(max_size(1)*mult));
  A.resize(floor(max_size(0)*mult), floor(max_size(1)*mult));

  MatrixXd Vf_scaled = d.V_uv_flat;

  for(unsigned i=0;i<3;++i)
    Vf_scaled.col(i) = Vf_scaled.col(i).array() - d.Vf_ori.col(i).minCoeff();
  for(unsigned i=0;i<2;++i)
    Vf_scaled.col(i) = Vf_scaled.col(i).array() / (d.Vf_ori.col(i).maxCoeff() - d.Vf_ori.col(i).minCoeff());
  for(unsigned i=0;i<2;++i)
    Vf_scaled.col(i) = (Vf_scaled.col(i).array() *2) - 1;


  Vf_scaled.col(2) = Vf_scaled.col(2).array() * 0;

  data.set_mesh(Vf_scaled, d.F_uv_flat);

  data.set_uv(d.V_uv, d.F_uv);


  data.set_colors(Eigen::RowVector3d(1,1,1));
  data.set_texture(d.R_final,d.G_final,d.B_final);

  core.proj = Matrix4f::Identity();
  core.model = Matrix4f::Identity();
  core.view = Matrix4f::Identity();

  core.background_color << 1.0f,1.0f,1.0f;

  // Do the rendering in a buffer
  core.draw_buffer(data,opengl,false,R,G,B,A);

  // Save the buffer in a PNG
  texture_to_png(img_filename, R, G, B);

}

void update_view(igl::viewer::Viewer& viewer, int mode)
{
  if (d.discoloration)
    d.update_discoloration_forward();

  if (mode == 1)
  {
    viewer.data.clear();
    viewer.data.set_mesh(d.V, d.F);

    viewer.core.show_texture = 0;

    if (d.V_uv.rows() > 0)
    {
      viewer.data.set_uv(d.V_uv, d.F_uv);
      viewer.core.show_texture = 1;
      viewer.data.set_texture(d.R,d.G,d.B);
    }

    viewer.data.set_colors(RowVector3d(1,1,1));

    viewer.core.show_lines = 0;

    MatrixXd C = MatrixXd::Zero(d.Vf.rows(),3);
    for (unsigned i=0; i<d.fixed.size(); ++i)
    {
      if (d.fixed(i) != 0)
        C.row(i) = RowVector3d(0,0,1);
      else
        C.row(i) = RowVector3d(1,0,0);
    }

    viewer.core.invert_normals = false;
    viewer.data.set_points(d.Vf, C);
  }

  if (mode == 2)
  {
    viewer.data.clear();
    MatrixXd Vtemp = MatrixXd::Zero(d.Vf.rows(),3);
    viewer.data.set_mesh(d.Vf, d.Ff);

    if (d.resample_uv)
      viewer.data.set_uv(d.Vf_uv_undistorted, d.Ff);
    else
    {
      MatrixXd tempuv = d.Vf_uv;
      for (unsigned i=0; i<2;++i)
        tempuv.col(i) = tempuv.col(i).array() - tempuv.col(i).minCoeff();

      // Scale to the correct size
      tempuv.col(0) = tempuv.col(0).array() / d.scale_u;
      tempuv.col(1) = tempuv.col(1).array() / d.scale_v;

      viewer.data.set_uv(tempuv, d.Ff);
    }

    viewer.data.set_texture(d.Rf,d.Gf,d.Bf);
    viewer.data.set_colors(Eigen::RowVector3d(1,1,1));
    viewer.core.show_texture = 1;

    viewer.core.invert_normals = true;

    MatrixXd Cpoint(d.Vf.rows(),3);
    for (unsigned i=0; i<Cpoint.rows(); ++i)
      Cpoint.row(i) = d.fixed(i) != 0 ? RowVector3d(1,0,0) : RowVector3d(0,1,0);

    viewer.data.set_points(d.Vf, Cpoint);
  }

  if (mode == 3)
  {
    viewer.data.clear();
    MatrixXd Vtemp = MatrixXd::Zero(4,3);
    Vtemp  << 0,0,0,
              1,0,0,
              1,1,0,
              0,1,0;
    MatrixXi Ftemp = MatrixXi::Zero(2,3);
    Ftemp << 0,1,2,2,3,0;

    viewer.data.set_uv(Vtemp, Ftemp);

    cerr << Vtemp << endl;

    Vector3d V_min = d.V.colwise().minCoeff();
    Vector3d V_max = d.V.colwise().maxCoeff();

    for (unsigned i=0;i<2;++i)
      Vtemp.col(i) = Vtemp.col(i).array() * (V_max(i)-V_min(i)) + V_min(i);

    cerr << Vtemp << endl;

    cerr << d.V.colwise().maxCoeff() << " " << d.V.colwise().minCoeff() << endl;

    viewer.data.set_mesh(Vtemp, Ftemp);
//    viewer.data.set_texture(d.in_triangulation.R,d.in_triangulation.G,d.in_triangulation.B);
    viewer.data.set_colors(Eigen::RowVector3d(1,1,1));
    viewer.core.show_texture = 1;
    viewer.core.show_lines = 0;
    viewer.core.align_camera_center(d.V,d.F);

    viewer.data.set_points(d.V, Eigen::RowVector3d(0,1,1));
  }

  if (mode == 4)
  {
    viewer.data.clear();

    MatrixXi Ftemp = d.Ff;

    for(unsigned i=0;i<Ftemp.rows();++i)
      if ((d.Vf(d.Ff(i,0),2) + d.Vf(d.Ff(i,1),2) + d.Vf(d.Ff(i,2),2)) > -0.00001)
        Ftemp.row(i) << 0,0,0;


    viewer.data.set_mesh(d.Vf, Ftemp);

    if (d.resample_uv)
      viewer.data.set_uv(d.Vf_uv_undistorted, d.Ff);
    else
    {
      MatrixXd tempuv = d.Vf_uv;
      for (unsigned i=0; i<2;++i)
        tempuv.col(i) = tempuv.col(i).array() - tempuv.col(i).minCoeff();

      // Scale to the correct size
      tempuv.col(0) = tempuv.col(0).array() / d.scale_u;
      tempuv.col(1) = tempuv.col(1).array() / d.scale_v;

      viewer.data.set_uv(tempuv, d.Ff);
    }

    viewer.data.set_texture(d.Rf,d.Gf,d.Bf);
    viewer.data.set_colors(Eigen::RowVector3d(1,1,1));
    viewer.core.show_texture = 1;
    viewer.core.show_lines = 0;
    viewer.core.invert_normals = true;

    viewer.data.set_points(d.V, Eigen::RowVector3d(1,0,0));
  }

  if (mode == 5)
  {
    d.invert_deformation();
    d.update_discoloration_inverse();

    viewer.data.clear();
    viewer.data.set_mesh(d.V_uv_flat, d.F_uv_flat);
    viewer.data.set_uv(d.V_uv, d.F_uv);

    viewer.data.set_texture(d.R_final,d.G_final,d.B_final);
    viewer.data.set_colors(Eigen::RowVector3d(1,1,1));
    viewer.core.show_texture = 1;
    viewer.core.show_lines = 0;
    viewer.core.invert_normals = true;

    if (d.Vq.rows() >0)
    {
      viewer.data.set_points(d.Vq_def, Eigen::RowVector3d(1,0,0));
      cerr << "Error: " << d.measure_error() << " mm" << "(" << d.measure_error()/(d.Fq.rows()*6) << ")" << endl;

    }

    save_image(std::string("temp_dipping.png"));
  }

  if (mode ==6)
  {
    viewer.data.clear();
    MatrixXd Vtemp = MatrixXd::Zero(d.Vf.rows(),3);
    viewer.data.set_mesh(d.Vf_uv, d.Ff);

    if (d.resample_uv)
      viewer.data.set_uv(d.Vf_uv_undistorted, d.Ff);
    else
    {
      MatrixXd tempuv = d.Vf_uv;
      for (unsigned i=0; i<2;++i)
        tempuv.col(i) = tempuv.col(i).array() - tempuv.col(i).minCoeff();

      // Scale to the correct size
      tempuv.col(0) = tempuv.col(0).array() / d.scale_u;
      tempuv.col(1) = tempuv.col(1).array() / d.scale_v;

      viewer.data.set_uv(tempuv, d.Ff);
    }

    viewer.data.set_texture(d.Rf,d.Gf,d.Bf);
    viewer.core.show_texture = 1;

  }

  if (mode ==7)
  {
    viewer.data.clear();
    viewer.data.set_mesh(d.V, d.F);

    viewer.core.show_texture = 0;
    viewer.data.set_colors(RowVector3d(0,0,0));
    viewer.core.show_lines = 0;
    viewer.core.invert_normals = false;
  }

  if (mode == 8)
  {
    viewer.data.clear();
    viewer.data.set_mesh(d.V, d.F);

    viewer.core.show_texture = 0;

    viewer.data.set_colors(RowVector3d(1,1,1));
    viewer.core.show_lines = 0;

    viewer.core.invert_normals = false;
  }

}

void save_cb(void *clientData)
{
  std::string fname = igl::file_dialog_save();
  if (fname.length() == 0)
    return;

  dipping_serialize(fname.c_str(), d);

  return;
}

void load_cb(void *clientData)
{
  std::string fname = igl::file_dialog_open();
  if (fname.length() == 0)
    return;

  dipping_deserialize(fname.c_str(), d);

  return;
}

void gridsearch_z_cb(void *clientData)
{
  int it = 10;
  for(unsigned i=0;i<it+1;++i)
  {
    double minv = 0.5;
    double maxv = 1;

    double v = minv + i*(maxv-minv)/double(it);

    d.stretch_multiplier = v;
    d.init_simulation();

    for (unsigned i=d.current_step; i<d.steps;++i)
    {
      cerr << i << endl;
      d.simulate();
    }

    d.invert_deformation();

    std::cerr.setf( std::ios::fixed, std:: ios::floatfield );
    std::cerr.precision(5);
    cerr << d.stretch_multiplier << "   " << ": " << d.measure_error() << " mm" << "(" << d.measure_error()/(d.Fq.rows()*6) << ")" << endl;

  }

}

IGL_INLINE void gridsearch_poisson_cb(void *clientData)
{
  int it = 10;
  for(unsigned i=0;i<it+1;++i)
  {
    double minv = 0;
    double maxv = 0.5;

    double v = minv + i*(maxv-minv)/double(it);

    d.poisson_ratio = v;
    d.init_simulation();

    for (unsigned i=d.current_step; i<d.steps;++i)
    {
      cerr << i << endl;
      d.simulate();
    }

    d.invert_deformation();

    std::cerr.setf( std::ios::fixed, std:: ios::floatfield );
    std::cerr.precision(5);
    cerr << d.poisson_ratio << "   " << ": " << d.measure_error() << " mm" << "(" << d.measure_error()/(d.Fq.rows()*6) << ")" << endl;

  }

}


IGL_INLINE void gridsearch_all_cb(void *clientData)
{

  char command[2048];
  sprintf(command, "mkdir -p %s", line_search_dir.c_str());
  system(command);


  VectorXd stretch_values, poisson_values, young_ratio_values;

  stretch_values.resize(11);
  for(unsigned i=0;i<11;++i)
    stretch_values[i] = 0.5+ i*(1.0-0.5)/10.;

  poisson_values.resize(2);
  poisson_values<<0., 0.001;//, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4;

  young_ratio_values.resize(9);
  young_ratio_values<< 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200, 500.0, 1000.0;

  int num = 0;
  char filename[2048];

  for(unsigned i=0;i<stretch_values.size();++i)
  {
    d.stretch_multiplier = stretch_values[i];

    for(unsigned j=0;j<poisson_values.size();++j)
    {
      d.poisson_ratio = poisson_values[j];

      for(unsigned k=0;k<young_ratio_values.size();++k, ++num)
      {
        d.min_young = 0.01;
        d.max_young = young_ratio_values[k] * d.min_young;
        d.do_adaptive_young = true;
        d.steps = 25;
        d.do_ARAP = false;
        d.resample_uv = true;
        d.init_simulation();

        printf("----- stretch: %.4g, poisson: %.4g, young max: %.4g ----- \n", d.stretch_multiplier, d.poisson_ratio, d.max_young);
        sprintf(filename, "%s/exp%d_0.txt", line_search_dir.c_str(),num);
        FILE *testfile = fopen(filename,"r");
        if (testfile)
        {
          cerr<<"Experiment "<<num<<" done."<<endl;
          fclose(testfile);
          continue;
        }

        bool ok = true;
        for (unsigned i=d.current_step; i<d.steps;++i)
        {
          if (!(d.simulate()))
          {
            ok = false;
            break;
          };
        }
        if (!ok)
        {
          cerr<<"error"<<endl;
          continue;
        }
        d.invert_deformation();

        std::cerr.setf( std::ios::fixed, std:: ios::floatfield );
        std::cerr.precision(5);
        double error = d.measure_error();
        cerr << error << " mm" << "(" << error/(d.Fq.rows()*6) << ")" << endl;


        ofstream ofs;
        sprintf(filename, "%s/exp%d_0.txt", line_search_dir.c_str(),num);
        ofs.open(filename);
        ofs<<stretch_values[i]<<" "<<poisson_values[j]<<" "<<young_ratio_values[k]<<" "<<error<<" "<<error/(d.Fq.rows()*6)<<endl;
        ofs.close();

      }
    }


  }

}


void run_simulation_for_new_model()
{
  char command[2048];
  sprintf(command, "mkdir -p %s", line_search_dir.c_str());
  system(command);

  //set IDf to something constant of the size of Rf, for the first iteration only
  d.IDf.resize(d.Rf.rows(),d.Rf.cols());
  d.IDf.setConstant((unsigned char)255);

  d.steps = 100;
  d.do_ARAP = false;
  d.resample_uv = true;
  d.stretch_multiplier = 0.65;
  d.do_adaptive_young = true;
  d.border_fixed = 1;

  d.zprojection = false;
  d.min_young = 0.01;
  d.max_young = 0.01;
  d.poisson_ratio = 0;
  double young_ratio = 600;

  bool converged = false;
  bool ok = true;
  char filename[1000];

  int num = 0;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Rc;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Gc;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Bc;

  while (!converged)
  {
    MatrixXd old_VFuv;
    if (d.resample_uv)
      old_VFuv = d.Vf;
    else
      old_VFuv = d.Vf_uv;
    //this computes the young values depending on IDf
    d.init_simulation();

    ok = true;
    for (unsigned i=d.current_step; i<d.steps;++i)
    {
      if (!(d.simulate()))
      {
        ok = false;
        break;
      };
    }
    if (!ok)
    {
      cerr<<"error"<<endl;
      break;
    }
    d.invert_deformation();
    d.update_discoloration_inverse();

    //save the image to be printed (on top of white bg for now)
    sprintf(filename, "%s/img%d.png",line_search_dir.c_str(), num);
    save_image(std::string(filename));

    texture_from_png(filename, R, G, B);
    d.Rf_ori = R;
    d.Gf_ori = G;
    d.Bf_ori = B;

    //set the IDf from d.R etc.
    d.set_ink_density_from_rgb(R, G, B);

    // Save the buffer in a PNG
    sprintf(filename, "%s/id%d.png",line_search_dir.c_str(), num);
    texture_to_png(filename, d.IDf, d.IDf, d.IDf);


    //set background to something with the desired density
    double avg_r=0., avg_g=0., avg_b=0.;
    int np = 0;
    for (unsigned i=0;i<R.rows();++i)
      for (unsigned j=0;j<R.cols();++j)
      {
        bool is_white = (R(i,j)==255 && G(i,j)==255 && B(i,j)==255);
        if (is_white)
          continue;
        avg_r += double(R(i,j))/255.0;
        avg_g += double(G(i,j))/255.0;
        avg_b += double(B(i,j))/255.0;
        np++;
      }
    avg_r/=np;
    avg_g/=np;
    avg_b/=np;
    double avg_c, avg_m, avg_y, avg_k;
    rgb2cmyk(avg_r, avg_g, avg_b, avg_c, avg_m, avg_y, avg_k);
    double bg_id =   0.2165;
    double s = avg_c+avg_m+avg_y+avg_k;
    avg_c *= bg_id/s;
    avg_m *= bg_id/s;
    avg_y *= bg_id/s;
    avg_k *= bg_id/s;
    cmyk2rgb(avg_c, avg_m, avg_y, avg_k, avg_r, avg_g, avg_b);
    Rc = R; Gc = G; Bc = B;
    for (unsigned i=0;i<R.rows();++i)
      for (unsigned j=0;j<R.cols();++j)
      {
        bool is_white = (R(i,j)==255 && G(i,j)==255 && B(i,j)==255);
        if (is_white)
        {
          Rc(i,j) = 255*avg_r;
          Gc(i,j) = 255*avg_g;
          Bc(i,j) = 255*avg_b;
        }
      }

    sprintf(filename, "%s/img_corr%d.png",line_search_dir.c_str(), num);
    texture_to_png(filename, Rc, Gc, Bc);



    d.max_young = young_ratio * d.min_young;
    if (num >0)
    {
      double diff;
      if (d.resample_uv)
        diff = (old_VFuv - d.Vf).rowwise().norm().mean();
      else
        diff = (old_VFuv - d.Vf_uv).rowwise().norm().mean();

      diff = diff / old_VFuv.rowwise().norm().mean();

      cerr<<endl<<endl<<"num = "<<num<<", diff = "<<diff<<endl;
      cerr<<"***********************************"<<endl<<endl;
      converged = diff<1e-4 || num >=5;
    }
    num++;
  }
  update_view(viewer,d.current_view);
}

bool callback_mouse_down(igl::viewer::Viewer& viewer, int button, int modifier)
{
  return false;
}

bool callback_key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
  if (key == '1')
  {
    d.current_view = 1;
    update_view(viewer,d.current_view);
  }

  if (key == '2')
  {
    d.current_view = 2;
    update_view(viewer,d.current_view);
  }

  if (key == '3')
  {
    d.current_view = 3;
    update_view(viewer,d.current_view);
  }

  if (key == '4')
  {
    d.current_view = 4;
    update_view(viewer,d.current_view);
  }

  if (key == '5')
  {
    d.current_view = 5;
    update_view(viewer,d.current_view);
  }

  if (key == '6')
  {
    d.current_view = 6;
    update_view(viewer,d.current_view);
  }

  if (key == '7')
  {
    d.current_view = 7;
    update_view(viewer,d.current_view);
  }

  if (key == '8')
  {
    d.current_view = 8;
    update_view(viewer,d.current_view);
  }

  if (key == '9')
  {
    d.current_view = 8;
    update_view(viewer,d.current_view);
  }

  if (key == '0')
  {
    d.current_view = 8;
    update_view(viewer,d.current_view);
  }

  if (key == 'Q')
  {
    d.init_simulation();
    update_view(viewer,d.current_view);
  }

  if (key == 'E')
  {
    d.simulate();
    update_view(viewer,d.current_view);
  }

  if (key == 'R')
  {
#ifdef TIME
    timer.start();
#endif

    for (unsigned i=d.current_step; i<d.steps;++i)
    {
      cerr << i << endl;
      d.simulate();
    }

    update_view(viewer,d.current_view);
//
//    if(d.discoloration)
//    {
//      d.update_discoloration_inverse_to_film();
//      d.Rf_ori = d.Rf;
//      d.Gf_ori = d.Gf;
//      d.Bf_ori = d.Bf;
//    }


#ifdef TIME
    timer.stop();
    cerr<<"total time (sec) "<<timer.getElapsedTimeInSec()<<endl;
#endif

  }

  if (key == 'N')
  {
    run_simulation_for_new_model();
  }

  if (key == 'T')
  {
    for (unsigned i=0; i<100;++i)
    {
      cerr << i << endl;
      d.simulate();
    }
    update_view(viewer,d.current_view);
  }

  if (key == 'O')
  {
    dipping_serialize("temp.libigl", d);
  }

  if (key == 'P')
  {
    dipping_deserialize("temp.libigl", d);
    update_view(viewer,d.current_view);
  }

  if (key == 'A')
  {
    d.discoloration = 0;
    char command[2048];
    sprintf(command, "mkdir -p %s", line_search_dir.c_str());
    system(command);
    for (unsigned i=d.current_step; i<d.steps;++i)
    {
      cerr << i << endl;
      d.simulate();

      char filename[1000];
      viewer.data.set_face_based(false);

      float lf;

      update_view(viewer,2);
      sprintf(filename, "%s/film_front%d.png", line_search_dir.c_str(),i);
      viewer.draw();
      igl::png::render_to_png(std::string(filename), viewer.core.viewport[2], viewer.core.viewport[3]);

    }

  }

  return true;
}

bool callback_init(igl::viewer::Viewer& viewer)
{
  d.init();
  d.init_simulation();

  // Init the new bar
  viewer.ngui->addWindow(Eigen::Vector2i(220,10),"Dipping");

  // add new group
  viewer.ngui->addGroup("IO");

  viewer.ngui->addButton("Load",[](){ load_cb(0); });
  viewer.ngui->addButton("Save",[](){ save_cb(0); });

  viewer.ngui->addGroup("Simulation");
  viewer.ngui->addButton("Init",[&viewer](){ callback_key_down(viewer, 'Q', 0);; });
  viewer.ngui->addButton("Step",[&viewer](){ callback_key_down(viewer, 'E', 0);; });
  viewer.ngui->addButton("All",[&viewer]() { callback_key_down(viewer, 'R', 0);; });

  viewer.ngui->addGroup("View");
  viewer.ngui->addButton("Obj",[&viewer](){ callback_key_down(viewer, '1', 0);; });
  viewer.ngui->addButton("Film",[&viewer](){ callback_key_down(viewer, '2', 0);; });
  viewer.ngui->addButton("Print",[&viewer]() { callback_key_down(viewer, '5', 0);; });

  viewer.ngui->addGroup("Params");
  viewer.ngui->addVariable("Steps",d.steps);
  viewer.ngui->addVariable("Stretch multiplier",d.stretch_multiplier);
  viewer.ngui->addVariable("Fix Border",d.border_fixed);
  viewer.ngui->addVariable("Resample_uv",d.resample_uv);
  viewer.ngui->addVariable("Z Projection",d.zprojection);
  viewer.ngui->addVariable("Discoloration",d.discoloration);
  viewer.ngui->addVariable("Discoloration Multiplier",d.discoloration_multiplier);

  viewer.ngui->addGroup("FEM");
  viewer.ngui->addVariable("do_ARAP",d.do_ARAP);
  viewer.ngui->addVariable("Inner iterations",d.inner_iter);
  viewer.ngui->addVariable("Young_modulus_nominal",d.young_modulus_nominal);
  viewer.ngui->addVariable("max_young",d.max_young);
  viewer.ngui->addVariable("min_young",d.min_young);
  viewer.ngui->addVariable("Adaptive Young",d.do_adaptive_young);
  viewer.ngui->addVariable("Poisson Ratio",d.poisson_ratio);
  viewer.ngui->addVariable("Symmetric FEM",d.fem_symmetric);
  viewer.ngui->addVariable("Plasticity",d.plasticity);

  viewer.ngui->addGroup("Grid Search");
  viewer.ngui->addButton("Z",[](){ gridsearch_z_cb(0); });
  viewer.ngui->addButton("Poisson",[](){ gridsearch_poisson_cb(0); });
  viewer.ngui->addButton("All",[](){ gridsearch_all_cb(0); });

  viewer.ngui->addGroup("Lighting");
  viewer.ngui->addVariable("Lighting",viewer.core.lighting_factor);

  viewer.screen->performLayout();
  return false;
}



template <typename DerivedV, typename DerivedF, typename DerivedT, typename Index>
IGL_INLINE bool readOBJPoly(
        const std::string str,
        Eigen::PlainObjectBase<DerivedV>& V,
        std::vector<std::vector< Index > >& F,
        Eigen::PlainObjectBase<DerivedV>& CN,
        Eigen::PlainObjectBase<DerivedF>& FN,
        Eigen::PlainObjectBase<DerivedT>& TC,
        Eigen::PlainObjectBase<DerivedF>& FTC)
{
  std::vector<std::vector<double> > vV,vTC,vN;
  std::vector<std::vector<Index> > vF,vFTC,vFN;
  bool success = igl::readOBJ(str,vV,vTC,vN,vF,vFTC,vFN);
  if(!success)
    return false;

  bool V_rect = igl::list_to_matrix(vV,V);
  if(!V_rect)
    return false;

  F = vF;

  if(!vN.empty())
  {
    bool VN_rect = igl::list_to_matrix(vN,CN);
    if(!VN_rect)
      return false;
  }

  if(!vFN.empty())
  {
    bool FN_rect = igl::list_to_matrix(vFN,FN);
    if(!FN_rect)
      return false;
  }

  if(!vTC.empty())
  {

    bool T_rect = igl::list_to_matrix(vTC,TC);
    if(!T_rect)
      return false;
  }
  if(!vFTC.empty())
  {

    bool FTC_rect = igl::list_to_matrix(vFTC,FTC);
    if(!FTC_rect)
      return false;
  }

  return true;
}

int main(int argc, char *argv[])
{
  if (!(argc == 4 || argc == 2|| argc == 5))
  {
    cout << "Usage iq_sim mesh.obj mesh.png plane.png" << endl;
    exit(0);
  }

  if (argc == 2)
  {
    dipping_deserialize(argv[1], d);
  }
  else
  {
    // Read scanned mesh ...
    Eigen::MatrixXd V,CN,TC;
    Eigen::MatrixXi F,FN,FTC;
    igl::readOBJ(argv[1],V,TC,CN,F,FTC,FN);

    d.V = V;
    d.F = F;

    d.V_uv = TC;
    d.F_uv = FTC;

    cerr << "Loading: " << argv[2] << endl;
    // ... and its texture
    texture_from_png(argv[2], d.R, d.G, d.B);


    // ... and the reference texture
    cerr << "Loading: " << argv[3] << endl;
    texture_from_png(argv[3], d.Rf_ori, d.Gf_ori, d.Bf_ori);
    texture_from_png(argv[3], d.Rf, d.Gf, d.Bf);

    // try to load the ink density image
    if (!texture_from_png(string(argv[3]) + "_ID.png", d.IDf, d.IDf, d.IDf))
    {
      cerr<<"Density image not found. Computing density..."<<endl;
      d.set_ink_density_from_rgb(d.Rf_ori, d.Gf_ori, d.Bf_ori);
    }
    if (argc>4)
      line_search_dir = std::string(argv[4]);


    if (d.IDf.rows() ==0)
      d.do_adaptive_young = false;


    // Loading the quad mesh
    MatrixXd Vq;
    vector<vector<int> > Fq;

    if (readOBJPoly(string(argv[1]) + ".quads.obj",Vq,Fq,CN,FN,TC,FTC))
    {
      d.Vq = Vq;
      d.Fq.resize(Fq.size(),4);
      for (unsigned i=0; i<Fq.size();++i)
        for (unsigned j=0; j<4;++j)
          d.Fq(i,j) = Fq[i][j];
    }


  }

  // Plot the mesh
  viewer.callback_key_down = callback_key_down;
  viewer.callback_mouse_down = callback_mouse_down;
  callback_key_down(viewer, '1', 0);
  viewer.core.align_camera_center(d.Vf,d.Ff);
  viewer.callback_init = callback_init;
  viewer.core.show_overlay = false;
  viewer.core.background_color << 1.0f,1.0f,1.0f;
  viewer.core.lighting_factor = 0.5;


  viewer.launch();
}

#else

using namespace Eigen;
using namespace std;
int main(int argc, char *argv[])
{
  return 0;

}
#endif
