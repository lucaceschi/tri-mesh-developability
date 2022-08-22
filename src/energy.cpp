#include "energy.hpp"

#include "mesh_matrix.hpp"


double combinatorialEnergy(const Eigen::Ref<const Matrix3Xd>& N,
                           const Eigen::Ref<const MatrixXi>& S,
                           const Eigen::Ref<const ArrayXb>& B)
{
    double totEnergy = 0.0;

    for(size_t v = 0; v < S.rows(); v++)
        totEnergy += localCombinatorialEnergy(v, N, S, B);

    return totEnergy;
}


double localCombinatorialEnergy(size_t v,
                                const Eigen::Ref<const Matrix3Xd>& N,
                                const Eigen::Ref<const MatrixXi>& S,
                                const Eigen::Ref<const ArrayXb>& B,
                                StarPartitioning* partitioning)
{
    int starSize;
    double energy = -1.0;

    for(starSize = 2; starSize < S.cols(); starSize++)
        if(S(v, starSize) == std::numeric_limits<std::size_t>::max())
            break;

    if(partitioning != nullptr)
    {
        partitioning->starSize = starSize;
        partitioning->v = v;
    }

    if(starSize <= 3 || B(v))
        energy = 0.0;
    else
    {
        // consider all possible cardinalities of a region resulting from a partitioning of the star
        for(int rSize = 2; rSize <= (starSize - 2); rSize++)
        {
            // consider all possible regions having rSize faces and not crossing the vector boundaries.
            // rBegin is the index of the first face in the region
            for(int rBegin = 0; rBegin < (starSize - rSize); rBegin++)
            {                
                double currRegionEnergy =  regionNormalDeviation(N, S, v, starSize, rBegin, rSize);
                double otherRegionEnergy = regionNormalDeviation(N, S, v, starSize, rBegin+rSize, starSize-rSize);
                
                double currPartitioningEnergy = currRegionEnergy + otherRegionEnergy;

                if(energy < 0 || currPartitioningEnergy < energy)
                {
                    energy = currPartitioningEnergy;
                    if(partitioning != nullptr)
                    {
                        partitioning->rBegin = rBegin;
                        partitioning->rSize = rSize;
                    }
                }
            }
        }
    }

    return energy;
}


double regionNormalDeviation(const Eigen::Ref<const Matrix3Xd>& N,
                             const Eigen::Ref<const MatrixXi>& S,
                             size_t v,
                             int starSize,
                             int rBegin,
                             int rSize)
{
    Eigen::Vector3d normDiff;
    double regionNormalDev = 0.0;

    for(int i = rBegin; i < (rBegin + rSize - 1); i++)
        for(int j = i+1; j < (rBegin + rSize); j++)
        {
            normDiff = ( N.row(S(v, i % starSize)) - N.row(S(v, j % starSize)) );
            regionNormalDev += (normDiff.squaredNorm() / std::pow(rSize, 2));
        }

    return regionNormalDev;
}