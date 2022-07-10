#ifndef MESH_MATRIX_HPP
#define MESH_MATRIX_HPP

#include "mesh.hpp"


using Matrix3Xd = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using Matrix3Xi = Eigen::Matrix<size_t, Eigen::Dynamic, 3, Eigen::RowMajor>;
using MatrixXi  = Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ArrayXd   = Eigen::ArrayXd;
using ArrayXb   = Eigen::Array<bool, Eigen::Dynamic, 1, Eigen::ColMajor>;


template <typename DerivedV, typename DerivedF>
extern void getMeshVF(const MyMesh& mesh,
                      Eigen::MatrixBase<DerivedV> const& V_,
                      Eigen::MatrixBase<DerivedF> const& F_);

template <typename DerivedS>
extern void getMeshStars(MyMesh& mesh,
                         Eigen::MatrixBase<DerivedS> const& S_);

template <typename DerivedB>
extern void getMeshBorders(MyMesh& mesh,
                           Eigen::ArrayBase<DerivedB> const& B_);

template <typename DerivedN, typename DerivedA>
extern void computeNormals(const Eigen::Ref<const Matrix3Xd>& V,
                           const Eigen::Ref<const Matrix3Xi>& F,
                           Eigen::MatrixBase<DerivedN> const& N_,
                           Eigen::ArrayBase<DerivedA> const& A_);

#endif