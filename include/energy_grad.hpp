#ifndef ENERGY_GRAD_HPP
#define ENERGY_GRAD_HPP

#include "mesh.hpp"
#include "energy.hpp"


/*
 * Compute the combinatorial energy of a mesh and its gradient
 * as defined in section B.2 of the paper "Developability of Triangle Meshes" by Stein et al. 2018
 */
double combinatorialEnergyGrad(MyMesh& m,
                             AreaFaceAttrHandle& fAttrArea,
                             StarVertAttrHandle& vAttrStar,
                             GradientVertAttrHandle& vAttrGrad);


/*
 * Compute the gradient of the normal deviation of a region
 * as defined in section B.2 of the paper "Developability of Triangle Meshes" by Stein et al. 2018
 */
void regionNormalDeviationGrad(MyMesh::VertexPointer v,
                               StarPartitioning& partitioning,
                               bool region,
                               MyMesh& m,
                               AreaFaceAttrHandle& fAttrArea,
                               StarVertAttrHandle& vAttrStar,
                               GradientVertAttrHandle& vAttrGrad);


/*
 * Compute the gradient of the normal of a face wrt one of its vertices
 * as defined in section B.1 of the paper "Developability of Triangle Meshes" by Stein et al. 2018
 */
vcg::Matrix33d faceNormalGrad(MyMesh::FacePointer f,
                              int vIndex,
                              AreaFaceAttrHandle& fAttrArea);


#endif