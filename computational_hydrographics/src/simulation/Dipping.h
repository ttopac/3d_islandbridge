#ifndef CH_DIPPING
#define CH_DIPPING

#include <Eigen/Core>
#include <vector>
#include <string>


using namespace std;
using namespace Eigen;

class Dipping
{
public:
  Dipping();
  ~Dipping();

  // Input mesh
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd V_uv;
  Eigen::MatrixXi F_uv;
  Eigen::MatrixXd V_uv_flat;
  Eigen::MatrixXi F_uv_flat;
  Eigen::VectorXi F_uv_flat_valid; //bool
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B;

  // Color enhanced texture for printing
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R_final;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G_final;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B_final;


  // Discolored Film
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Rd;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Gd;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Bd;

  // Texture mesh
  Eigen::MatrixXd Vf;
  Eigen::MatrixXd Vf_ori;
  Eigen::MatrixXi Ff;
  Eigen::MatrixXd Vf_uv;

  //Strain values
  Eigen::Matrix2d S;
  Eigen::Matrix2d I;

  //Olga: are undistorted and ori the same?
  Eigen::MatrixXd Vf_uv_undistorted;
  Eigen::MatrixXd Vf_uv_ori;
  Eigen::VectorXi Vf_uv_fid;
  Eigen::MatrixXi TTf;
  Eigen::MatrixXi TTif;
  double scale_u;
  double scale_v;

  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Rf_ori;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Gf_ori;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Bf_ori;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Rf;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Gf;
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> Bf;

  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IDf;//ink density

  Eigen::VectorXi fixed;
  Eigen::VectorXi new_fixed;
  vector<bool> border;

  // Quad mesh for validation
  Eigen::MatrixXd Vq;
  Eigen::MatrixXd Vq_bc;
  Eigen::MatrixXd Vq_def;
  Eigen::MatrixXi Fq;

  // Simulation parameters
  double dt;
  int current_view;
  double steps;
  double stretch_multiplier;
  double current_step;
  bool border_fixed;
  bool do_ARAP;
  bool resample_uv;
  //todo:serialize, FEM stuff
  int inner_iter;
  double young_modulus_nominal;
  Eigen::VectorXd young_modulus_array;//adaptive young modulus per trianlge;
  bool do_adaptive_young;
  double min_young, max_young;
  //todo: binning int num_young_bins;

  double poisson_ratio;
  bool fem_symmetric;
  bool plasticity;

  bool zprojection;
  bool discoloration;
  double discoloration_multiplier;

  // Initialize the texture mesh
  void init_texture_mesh(double squaresize, int u_count, int v_count);

  // Calculate the per-face young modulus values, based on ink density of texture image
  void calculate_young_modulus_values();

  // Compute the ink density image from RGB image
  void set_ink_density_from_rgb(const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,
                                const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G,
                                const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B);

  // General initialization
  void init();

  // Initialize the simulation
  void init_simulation();

  // Main simulation function, true if simulation is finished
  bool simulate();

  // Resample the mesh using the UV coordinates
  void resample();

  // Substeps, do not call
  void move_object_and_film_down_project();
  bool deform();
  bool deformARAP();
  bool deformFEM();

  // project the scanned mesh into the plane
  void invert_deformation();

  // compute a quality measure of the current state
  double measure_error();

  // compute discoloration
  void update_discoloration_forward();
  void update_discoloration_inverse();
  void update_discoloration_inverse_to_film();

  // Debug
  MatrixXd C_vert;

  void fill_black(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B);
};

#endif
