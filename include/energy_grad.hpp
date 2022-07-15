#ifndef ENERGY_GRAD_HPP
#define ENERGY_GRAD_HPP

#include "mesh.hpp"
#include "energy.hpp"


// -------------------------------------------------------------------------------------------------
// Energy gradient on matrix representation

Eigen::Matrix3d faceNormalGrad(size_t f,
                               int FVindex,
                               const Eigen::Ref<const Matrix3Xd>& V,
                               const Eigen::Ref<const Matrix3Xi>& F,
                               const Eigen::Ref<const Matrix3Xd>& N,
                               const Eigen::Ref<const ArrayXd>& A);

void regionNormalDeviationGrad(StarPartitioning& region,
                               const Eigen::Ref<const Matrix3Xd>& V,
                               const Eigen::Ref<const Matrix3Xi>& F,
                               const Eigen::Ref<const Matrix3Xd>& N,
                               const Eigen::Ref<const ArrayXd>& A,
                               const Eigen::Ref<const MatrixXi>& S,
                               Eigen::Ref<Matrix3Xd> G);

void combinatorialEnergyGrad(const Eigen::Ref<const Matrix3Xd>& V,
                             const Eigen::Ref<const Matrix3Xi>& F,
                             const Eigen::Ref<const Matrix3Xd>& N,
                             const Eigen::Ref<const ArrayXd>& A,
                             const Eigen::Ref<const MatrixXi>& S,
                             const Eigen::Ref<const ArrayXb>& B,
                             Eigen::Ref<Matrix3Xd> G,
                             std::function<void(double, size_t)> localEnergyCallback);

void combinatorialEnergyGrad(const Eigen::Ref<const Matrix3Xd>& V,
                             const Eigen::Ref<const Matrix3Xi>& F,
                             const Eigen::Ref<const Matrix3Xd>& N,
                             const Eigen::Ref<const ArrayXd>& A,
                             const Eigen::Ref<const MatrixXi>& S,
                             const Eigen::Ref<const ArrayXb>& B,
                             Eigen::Ref<Matrix3Xd> G);

#endif