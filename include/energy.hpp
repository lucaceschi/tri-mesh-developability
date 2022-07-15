#ifndef ENERGY_HPP
#define ENERGY_HPP

#include <vector>

#include "mesh.hpp"
#include "mesh_matrix.hpp"


struct StarPartitioning
{
    size_t v;
    int starSize;
    int rBegin;   // index of the starting face within star for the first region
    int rSize;    // cardinality of the first region
};

double regionNormalDeviation(const Eigen::Ref<const Matrix3Xd>& N,
                             const Eigen::Ref<const MatrixXi>& S,
                             size_t v,
                             int starSize,
                             int rBegin,
                             int rSize);

double localCombinatorialEnergy(size_t v,
                                const Eigen::Ref<const Matrix3Xd>& N,
                                const Eigen::Ref<const MatrixXi>& S,
                                const Eigen::Ref<const ArrayXb>& B,
                                StarPartitioning* partitioning = nullptr);

double combinatorialEnergy(const Eigen::Ref<const Matrix3Xd>& N,
                           const Eigen::Ref<const MatrixXi>& S,
                           const Eigen::Ref<const ArrayXb>& B);



#endif