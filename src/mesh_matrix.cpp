#include "mesh_matrix.hpp"

#include <limits>


template <typename DerivedV, typename DerivedF>
void getMeshVF(const MyMesh& mesh,
                    Eigen::MatrixBase<DerivedV> const& V_,
                    Eigen::MatrixBase<DerivedF> const& F_)
{    
    Eigen::MatrixBase<DerivedV>& V = const_cast< Eigen::MatrixBase<DerivedV>& >(V_);
    Eigen::MatrixBase<DerivedF>& F = const_cast< Eigen::MatrixBase<DerivedF>& >(F_);
    
    vcg::tri::RequireCompactness(mesh);

    V.derived().resize(mesh.VN(), 3);
    for(size_t i = 0; i < mesh.VN(); i++)
        for(size_t j = 0; j < 3; j++)
            V(i,j) = mesh.vert[i].cP()[j];

    F.derived().resize(mesh.FN(), 3);
    for(size_t i = 0; i < mesh.FN(); i++)
        for(size_t j = 0; j < 3; j++)
            F(i,j) = (size_t) vcg::tri::Index(mesh, mesh.face[i].cV(j));
}

template void getMeshVF<Matrix3Xd, Matrix3Xi>(
    const MyMesh& mesh,
    Eigen::MatrixBase<Matrix3Xd> const& V_,
    Eigen::MatrixBase<Matrix3Xi> const& F_
);


template <typename DerivedS>
void getMeshStars(MyMesh& mesh,
                  Eigen::MatrixBase<DerivedS> const& S_)
{
    Eigen::MatrixBase<DerivedS>& S = const_cast< Eigen::MatrixBase<DerivedS>& >(S_);

    MyMesh::VertexPointer vPointer;
    std::vector<vcg::face::Pos<MyFace>> currStarPos;
    std::vector<int> starSizes;
    size_t vIndex;
    int i;

    S.derived().resize(mesh.VN(), 2);
    starSizes.reserve(mesh.VN());
    for(vIndex = 0; vIndex < mesh.VN(); vIndex++) {
        vPointer = &mesh.vert[vIndex];
        vcg::face::VFOrderedStarFF(vcg::face::Pos<MyFace>(vPointer->VFp(), vPointer), currStarPos);
        starSizes.push_back(currStarPos.size());

        if(currStarPos.size() > S.cols())
            S.derived().conservativeResize(Eigen::NoChange, currStarPos.size());

        for(i = 0; i < S.cols(); i++)
            S(vIndex, i) = vcg::tri::Index(mesh, currStarPos[i].F());
    }

    for(vIndex = 0; vIndex < mesh.VN(); vIndex++)
        for(i = starSizes[vIndex]; i < S.cols(); i++)
            S(vIndex, i) = std::numeric_limits<std::size_t>::max();
}

template void getMeshStars<MatrixXi>(
    MyMesh& mesh,
    Eigen::MatrixBase<MatrixXi> const& S_
);


template <typename DerivedB>
void getMeshBorders(MyMesh& mesh,
                    Eigen::ArrayBase<DerivedB> const& B_)
{
    Eigen::ArrayBase<DerivedB>&  B = const_cast< Eigen::ArrayBase<DerivedB>& >(B_);

    B.derived().resize(mesh.VN(), Eigen::NoChange);
    for(size_t vIndex = 0; vIndex < B.rows(); vIndex++)
        B(vIndex) = mesh.vert[vIndex].IsB();
}

template void getMeshBorders<ArrayXb>(
    MyMesh& mesh,
    Eigen::ArrayBase<ArrayXb> const& B_
);


template <typename DerivedN, typename DerivedA>
void computeNormals(const Eigen::Ref<const Matrix3Xd>& V,
                    const Eigen::Ref<const Matrix3Xi>& F,
                    Eigen::MatrixBase<DerivedN> const& N_,
                    Eigen::ArrayBase<DerivedA> const& A_)
{
    typedef typename DerivedA::Scalar AScalar;
    Eigen::MatrixBase<DerivedN>& N = const_cast< Eigen::MatrixBase<DerivedN>& >(N_);
    Eigen::ArrayBase<DerivedA>&  A = const_cast< Eigen::ArrayBase<DerivedA>&  >(A_);
    
    Eigen::RowVector3d edgeAvec, edgeBvec;
    AScalar norm;
    
    N.derived().resize(F.rows(), Eigen::NoChange);
    A.derived().resize(F.rows());
    for(size_t fIndex = 0; fIndex < F.rows(); fIndex++)
    {
        edgeAvec = V.row(F(fIndex, 1)) - V.row(F(fIndex, 0));
        edgeBvec = V.row(F(fIndex, 2)) - V.row(F(fIndex, 0));

        N.row(fIndex) = edgeAvec.cross(edgeBvec);
        norm = N.row(fIndex).norm();
        A(fIndex) = norm / 2.0;
        assert(norm > 0); // the two edges cannot be parallel
        N.row(fIndex) /= norm;
    }
}

template void computeNormals<Matrix3Xd, ArrayXd>(
    const Eigen::Ref<const Matrix3Xd>& V,
    const Eigen::Ref<const Matrix3Xi>& F,
    Eigen::MatrixBase<Matrix3Xd> const& N_,
    Eigen::ArrayBase<ArrayXd> const& A_
);