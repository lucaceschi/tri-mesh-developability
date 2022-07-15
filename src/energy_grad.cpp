#include "energy_grad.hpp"

#include "mesh_matrix.hpp"


Eigen::Matrix3d faceNormalGrad(size_t f,
                               int FVindex,
                               const Eigen::Ref<const Matrix3Xd>& V,
                               const Eigen::Ref<const Matrix3Xi>& F,
                               const Eigen::Ref<const Matrix3Xd>& N,
                               const Eigen::Ref<const ArrayXd>& A)
{
    Eigen::Vector3d oppositeEdge = V.row(F(f, (FVindex+2)%3)) - V.row(F(f, (FVindex+1)%3));
    Eigen::Vector3d normal = N.row(f);
    Eigen::Matrix3d grad = (oppositeEdge.cross(normal) * normal.transpose());

    return grad / A(f);
}


void regionNormalDeviationGrad(StarPartitioning& region,
                               const Eigen::Ref<const Matrix3Xd>& V,
                               const Eigen::Ref<const Matrix3Xi>& F,
                               const Eigen::Ref<const Matrix3Xd>& N,
                               const Eigen::Ref<const ArrayXd>& A,
                               const Eigen::Ref<const MatrixXi>& S,
                               Eigen::Ref<Matrix3Xd> G)
{  
    size_t faceA, faceB;
    Eigen::Vector3d normalDiff;

    int faceA_FVindex, faceB_FVindex;
    size_t vertA, vertB;
    bool sharedV;
    
    for(int i = region.rBegin; i < (region.rBegin + region.rSize - 1); i++)
        for(int j = i+1; j < (region.rBegin + region.rSize); j++)
        {
            faceA = S(region.v, i % region.starSize);
            faceB = S(region.v, j % region.starSize);

            normalDiff = N.row(faceA) - N.row(faceB);

            for(faceA_FVindex = 0; faceA_FVindex < 3; faceA_FVindex++)
            {
                vertA = F(faceA, faceA_FVindex);
                G.row(vertA) += (normalDiff.transpose() * faceNormalGrad(faceA, faceA_FVindex, V, F, N, A) * 2/std::pow(region.rSize, 2));
            }

            for(faceB_FVindex = 0; faceB_FVindex < 3; faceB_FVindex++)
            {
                vertB = F(faceB, faceB_FVindex);
                G.row(vertB) -= (normalDiff.transpose() * faceNormalGrad(faceB, faceB_FVindex, V, F, N, A) * 2/std::pow(region.rSize, 2));
            }            
        }
}


void combinatorialEnergyGrad(const Eigen::Ref<const Matrix3Xd>& V,
                             const Eigen::Ref<const Matrix3Xi>& F,
                             const Eigen::Ref<const Matrix3Xd>& N,
                             const Eigen::Ref<const ArrayXd>& A,
                             const Eigen::Ref<const MatrixXi>& S,
                             const Eigen::Ref<const ArrayXb>& B,
                             Eigen::Ref<Matrix3Xd> G,
                             std::function<void(double, size_t)> localEnergyCallback)
{
    G.setZero();

    double currEnergy;
    StarPartitioning currPart;
    for(size_t v = 0; v < V.rows(); v++)
    {
        currEnergy = localCombinatorialEnergy(v, N, S, B, &currPart);
        localEnergyCallback(currEnergy, v);

        if(currPart.starSize <= 3 || B(v))
            continue;

        regionNormalDeviationGrad(currPart, V, F, N, A, S, G);
        currPart.rBegin = (currPart.rBegin + currPart.rSize); // % currPart.starI.size();
        currPart.rSize  = currPart.starSize - currPart.rSize;
        regionNormalDeviationGrad(currPart, V, F, N, A, S, G);
    }
}


void combinatorialEnergyGrad(const Eigen::Ref<const Matrix3Xd>& V,
                             const Eigen::Ref<const Matrix3Xi>& F,
                             const Eigen::Ref<const Matrix3Xd>& N,
                             const Eigen::Ref<const ArrayXd>& A,
                             const Eigen::Ref<const MatrixXi>& S,
                             const Eigen::Ref<const ArrayXb>& B,
                             Eigen::Ref<Matrix3Xd> G)
{
    combinatorialEnergyGrad(V, F, N, A, S, B, G, [](double e, size_t v) {});
}