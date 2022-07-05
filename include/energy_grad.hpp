#ifndef ENERGY_GRAD_HPP
#define ENERGY_GRAD_HPP

#include "mesh.hpp"
#include "energy.hpp"


vcg::Matrix33d faceNormalGrad(MyFace* f, int v, MyMesh::PerFaceAttributeHandle<double>& triArea);

void combinatorialEnergyGrad(MyMesh& mesh,
                             MyMesh::PerVertexAttributeHandle<vcg::Point3d> &energyGrad,
                             std::function<void(double, MyVertex*)> localEnergyCallback);

void regionNormalDeviationGrad(MyVertex* v,
                               std::vector<vcg::face::Pos<MyFace>> star,
                               int regionBeginIndex,
                               int regionSize,
                               MyMesh::PerVertexAttributeHandle<vcg::Point3d>& energyGrad);

#endif