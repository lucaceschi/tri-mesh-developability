#ifndef MESH_MATRIX_HPP
#define MESH_MATRIX_HPP

#include "mesh.hpp"


using Matrix3Xd = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using Matrix3Xi = Eigen::Matrix<size_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using MatrixXi  = Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ArrayXd   = Eigen::ArrayXd;
using ArrayXb   = Eigen::Array<bool, Eigen::Dynamic, 1, Eigen::ColMajor>;


void getMeshVF(const MyMesh& mesh,
               Eigen::Ref<Matrix3Xd> V,
               Eigen::Ref<Matrix3Xi> F);

template <typename DerivedS>
extern void getMeshStars(MyMesh& mesh,
                         Eigen::MatrixBase<DerivedS> const& S_);

void getMeshBorders(MyMesh& mesh,
                    Eigen::Ref<ArrayXb> B);

void computeNormals(const Eigen::Ref<const Matrix3Xd>& V,
                    const Eigen::Ref<const Matrix3Xi>& F,
                    Eigen::Ref<Matrix3Xd> N,
                    Eigen::Ref<ArrayXd> A);

#endif