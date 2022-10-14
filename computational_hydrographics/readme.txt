This is the reference implementation of:

Texture Mapping Real-World Objects with Hydrographics
Daniele Panozzo, Olga Diamanti, Sylvain Paris, Marco Tarini, Evgeni Sorkine, Olga Sorkine-Hornung

as presented at SGP 2015. This code was used to create all the results in the paper, the only difference is that in this version we replaced the commercial solver PARDISO with the open source sparse LU solver in Eigen. This does not affect the results, but it considerably reduces performances on workstations
with multiple cores. The advantage is that all the source code in this package
is open source.

Licence: The code is distributed with the MPL2 licence. If you use it in a scientific publication or in a product, please cite our paper:

@article{Hydrographics:Panozzo:2015,
author = {Daniele Panozzo and Olga Diamanti and Sylvain Paris and Marco Tarini and Evgeni Sorkine and Olga Sorkine-Hornung},
title = {Texture Mapping Real-World Objects with Hydrographics},
journal = {Computer Graphics Forum (proceedings of EUROGRAPHICS Symposium on Geometry Processing)},
year = {2015},
}

IMPORTANT: The data in the folder "data" can only be used for comparison purposes in a scientific publication. You cannot use this data for any other purpose without a written authorization from Daniele Panozzo.

Compilation:

The code has only been tested on MacOSX 10.10 but it should be easy to port it to linux or windows. It depends on libigl, which should be downloaded in the same folder where the zip file containing this readme has been unpacked. The directory layout should be:

---- parent folder
       |----- libigl
       |----- computational_hydrographics

To download libigl:

git clone https://github.com/libigl/libigl.git
cd libigl
git submodule update --init --recursive
[Optional: git checkout 0022da9] <-- this is the version of libigl I tested it with.

After libigl is downloaded, compile embree:

cd libigl/external/embree
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make

(this step might fail if you do not have ispc installed, in that case download
the latest version from the ispc website and copy it in the build directory before
building embree)

Finally compile the hydrographics prototype:

cd computational_hydrographics
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make


How to use it:

./ig_sim data/Ac.obj data/A.png checkerboard_light.png

and similarly with the other provided datasets. Loading datasets might take up to a minute, wait for a while if you see a black screen after issuing this command before killing it.

Then press Simulation -> All. The result will be computed in a few minutes
and you can then observe it using the three visualization options:

View->Obj: shows the object being dipped
View->Film: shows the film durin the immersion
View->Print: shows the pattern to print to obtain the desidered texture

If you want to speed up the algorithm to see a preview of the result, change the number of iterations to 10, press Simulation->Init and then Simulation->All. You
can also see the intermediate steps by pressing Simulation->Init and then
multiple times Simulation->Step.

Copyright Daniele Panozzo (daniele.panozzo@gmail.com) and Olga Diamanti, 2015
