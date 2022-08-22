#ifndef ENERGY_GRAD_HPP
#define ENERGY_GRAD_HPP

#include "mesh.hpp"
#include "energy.hpp"


vcg::Matrix33d faceNormalGrad(MyMesh::FacePointer f,
                              int vIndex,
                              MyMesh& m,
                              AreaFaceAttrHandle& fAttrArea);

void regionNormalDeviationGrad(MyMesh::VertexPointer v,
                               StarPartitioning& partitioning,
                               bool region,
                               MyMesh& m,
                               AreaFaceAttrHandle& fAttrArea,
                               StarVertAttrHandle& vAttrStar,
                               GradientVertAttrHandle& vAttrGrad);

double combinatorialEnergyGrad(MyMesh& m,
                             AreaFaceAttrHandle& fAttrArea,
                             StarVertAttrHandle& vAttrStar,
                             GradientVertAttrHandle& vAttrGrad);

#endif