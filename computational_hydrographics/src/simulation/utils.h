//
//  utils.h
//  ig
//
//  Created by Olga Diamanti on 24/03/15.
//
//

#ifndef __ig__utils__
#define __ig__utils__

#include <Eigen/Core>
void rgb2cmyk(double r, double g, double b, double& c,double& m,double& y,double& k);
void cmyk2rgb(double c,double m,double y,double k, double& r, double& g, double& b);
bool texture_from_png(
                      const std::string png_file,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
                      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B
                      );
bool texture_to_png(
                    const std::string png_file,
                    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R_,
                    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G_,
                    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B_
                    );

#endif /* defined(__ig__utils__) */
