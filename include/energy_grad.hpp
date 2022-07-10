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

template <typename DerivedG>
extern void regionNormalDeviationGrad(StarPartitioning& region,
                                      const Eigen::Ref<const Matrix3Xd>& V,
                                      const Eigen::Ref<const Matrix3Xi>& F,
                                      const Eigen::Ref<const Matrix3Xd>& N,
                                      const Eigen::Ref<const ArrayXd>& A,
                                      const Eigen::Ref<const MatrixXi>& S,
                                      Eigen::MatrixBase<DerivedG> const& G_);

template <typename DerivedG>
extern void combinatorialEnergyGrad(const Eigen::Ref<const Matrix3Xd>& V,
                                    const Eigen::Ref<const Matrix3Xi>& F,
                                    const Eigen::Ref<const Matrix3Xd>& N,
                                    const Eigen::Ref<const ArrayXd>& A,
                                    const Eigen::Ref<const MatrixXi>& S,
                                    const Eigen::Ref<const ArrayXb>& B,
                                    Eigen::MatrixBase<DerivedG> const& G_,
                                    std::function<void(double, size_t)> localEnergyCallback);

#endif