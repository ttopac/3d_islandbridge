#include "Dipping_serialize.h"
#include <igl/is_border_vertex.h>

//#define SER1(X)  igl::serialize(d.X,"X",filename,true);
//#define SER(X)   igl::serialize(d.X,"X",filename,false);
//#define DESER(X) igl::deserialize(d.X,"X",filename);

void dipping_serialize(const std::string filename, const Dipping& d)
{
  igl::serialize(d.V,"V",filename,true);
  igl::serialize(d.F,"F",filename);

  igl::serialize(d.V_uv,"V_uv",filename);
  igl::serialize(d.F_uv,"F_uv",filename);

  igl::serialize(d.V_uv_flat,"V_uv_flat",filename);
  igl::serialize(d.F_uv_flat,"F_uv_flat",filename);
  igl::serialize(d.R,"R",filename);
  igl::serialize(d.G,"G",filename);
  igl::serialize(d.B,"B",filename);
  igl::serialize(d.R,"Rd",filename);
  igl::serialize(d.G,"Gd",filename);
  igl::serialize(d.B,"Bd",filename);
  igl::serialize(d.R,"R_final",filename);
  igl::serialize(d.G,"G_final",filename);
  igl::serialize(d.B,"B_final",filename);


  igl::serialize(d.Vf,"Vf",filename);
  igl::serialize(d.Vf_ori,"Vf_ori",filename);
  igl::serialize(d.Ff,"Ff",filename);
  igl::serialize(d.Vf_uv,"Vf_uv",filename);
  igl::serialize(d.Vf_uv_undistorted,"Vf_uv_undistorted",filename);
  igl::serialize(d.Vf_uv_ori,"Vf_uv_ori",filename);
  igl::serialize(d.Vf_uv_fid,"Vf_uv_fid",filename);
  igl::serialize(d.TTf,"TTf",filename);
  igl::serialize(d.TTif,"TTif",filename);

  igl::serialize(d.scale_u,"scale_u",filename);
  igl::serialize(d.scale_v,"scale_v",filename);
  igl::serialize(d.Rf,"Rf",filename);
  igl::serialize(d.Gf,"Gf",filename);
  igl::serialize(d.Bf,"Bf",filename);
  igl::serialize(d.fixed,"fixed",filename);

  igl::serialize(d.dt,"dt",filename);
  igl::serialize(d.current_view,"current_view",filename);
  igl::serialize(d.steps,"steps",filename);
  igl::serialize(d.stretch_multiplier,"stretch_multiplier",filename);
  igl::serialize(d.current_step,"current_step",filename);

  igl::serialize(d.border_fixed,"border_fixed",filename);
  igl::serialize(d.resample_uv,"resample_uv",filename);

  igl::serialize(d.Vq,"Vq",filename);
  igl::serialize(d.Fq,"Fq",filename);
  igl::serialize(d.Vq_bc,"Vq_bc",filename);
  igl::serialize(d.Vq_def,"Vq_def",filename);

  igl::serialize(d.border,"border",filename);
//  igl::serialize(d.FEM_W,"FEM_W",filename);
  
  igl::serialize(d.zprojection,"zprojection",filename);
  igl::serialize(d.discoloration,"discoloration",filename);
  igl::serialize(d.discoloration_multiplier,"discoloration_multiplier",filename);


}

void dipping_deserialize(const std::string filename, Dipping& d)
{
  igl::deserialize(d.V,"V",filename);
  igl::deserialize(d.F,"F",filename);

  igl::deserialize(d.V_uv,"V_uv",filename);
  igl::deserialize(d.F_uv,"F_uv",filename);

  igl::deserialize(d.V_uv_flat,"V_uv_flat",filename);
  igl::deserialize(d.F_uv_flat,"F_uv_flat",filename);
  igl::deserialize(d.R,"R",filename);
  igl::deserialize(d.G,"G",filename);
  igl::deserialize(d.B,"B",filename);
  igl::deserialize(d.R,"Rd",filename);
  igl::deserialize(d.G,"Gd",filename);
  igl::deserialize(d.B,"Bd",filename);
  igl::deserialize(d.R,"R_final",filename);
  igl::deserialize(d.G,"G_final",filename);
  igl::deserialize(d.B,"B_final",filename);

  igl::deserialize(d.Vf,"Vf",filename);
  igl::deserialize(d.Vf_ori,"Vf_ori",filename);
  igl::deserialize(d.Ff,"Ff",filename);
  igl::deserialize(d.Vf_uv,"Vf_uv",filename);
  igl::deserialize(d.Vf_uv_undistorted,"Vf_uv_undistorted",filename);
  igl::deserialize(d.Vf_uv_ori,"Vf_uv_ori",filename);
  igl::deserialize(d.Vf_uv_fid,"Vf_uv_fid",filename);
  igl::deserialize(d.TTf,"TTf",filename);
  igl::deserialize(d.TTif,"TTif",filename);

  igl::deserialize(d.scale_u,"scale_u",filename);
  igl::deserialize(d.scale_v,"scale_v",filename);
  igl::deserialize(d.Rf,"Rf",filename);
  igl::deserialize(d.Gf,"Gf",filename);
  igl::deserialize(d.Bf,"Bf",filename);
  igl::deserialize(d.fixed,"fixed",filename);

  igl::deserialize(d.dt,"dt",filename);
  igl::deserialize(d.current_view,"current_view",filename);
  igl::deserialize(d.steps,"steps",filename);
  igl::deserialize(d.stretch_multiplier,"stretch_multiplier",filename);
  igl::deserialize(d.current_step,"current_step",filename);

  igl::deserialize(d.border_fixed,"border_fixed",filename);
  igl::deserialize(d.resample_uv,"resample_uv",filename);

  igl::deserialize(d.Vq,"Vq",filename);
  igl::deserialize(d.Fq,"Fq",filename);
  igl::deserialize(d.Vq_bc,"Vq_bc",filename);
  igl::deserialize(d.Vq_def,"Vq_def",filename);
  
  igl::deserialize(d.border,"border",filename);
  
//  igl::deserialize(d.FEM_W,"FEM_W",filename);

  igl::deserialize(d.zprojection,"zprojection",filename);
  igl::deserialize(d.discoloration,"discoloration",filename);
  igl::deserialize(d.discoloration_multiplier,"discoloration_multiplier",filename);

  //d.border = igl::is_border_vertex(d.Vf,d.Ff);


}
