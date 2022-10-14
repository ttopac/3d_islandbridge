//
//  utils.cpp
//  ig
//
//  Created by Olga Diamanti on 24/03/15.
//
//

#include "utils.h"
#include "../../ext/png/png.hpp"
#include <algorithm>
void rgb2cmyk(double r, double g, double b, double& c,double& m,double& y,double& k)
{
  k = 1.-std::max(std::max(r,g),b);
  if (k >= 1)
  {
    c = 0;
    m = 0;
    y = 0;
  }
  else
  {
    c = (1.-r-k) / (1.-k);
    m = (1.-g-k) / (1.-k);
    y = (1.-b-k) / (1.-k);
  }
}

void cmyk2rgb(double c,double m,double y,double k, double& r, double& g, double& b)
{
  r = (1.-c) * (1.-k);
  g = (1.-m) * (1.-k);
  b = (1.-y) * (1.-k);
}

bool texture_from_png(
                      const std::string png_file,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B
                      )
{
  FILE *fp = fopen(png_file.c_str(), "r");
  if(!fp)
    return false;
  else
    fclose(fp);
    
  png::image< png::rgb_pixel > image(png_file);
  
  R.resize(image.get_height(),image.get_width());
  G.resize(image.get_height(),image.get_width());
  B.resize(image.get_height(),image.get_width());
  
  unsigned h = image.get_height();
  unsigned w = image.get_width();
  
  for (unsigned i=0; i<h; ++i)
  {
    for (unsigned j=0; j<w; ++j)
    {
      R(i,j) = image.get_pixel(j,h-1-i).red;
      G(i,j) = image.get_pixel(j,h-1-i).green;
      B(i,j) = image.get_pixel(j,h-1-i).blue;
    }
  }
  
  R.transposeInPlace();
  G.transposeInPlace();
  B.transposeInPlace();
  
  return true;
}

bool texture_to_png(
                    const std::string png_file,
                    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R_,
                    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G_,
                    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B_
                    )
{
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R = R_.transpose();
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G = G_.transpose();
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B = B_.transpose();
  
  unsigned h = R.rows();
  unsigned w = R.cols();
  
  png::image< png::rgb_pixel > image(w,h);
  
  for (unsigned i=0; i<h; ++i)
  {
    for (unsigned j=0; j<w; ++j)
    {
      image.set_pixel(j,h-1-i,png::rgb_pixel(R(i,j), G(i,j), B(i,j)));
    }
  }
  
  image.write(png_file);
  
  return true;
}