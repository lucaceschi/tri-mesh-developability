#include "energy.hpp"


double regionNormalDeviation(std::vector<vcg::face::Pos<MyFace>> star,
                             int regionBeginIndex,
                             int regionSize)
{   
    vcg::Point3d n;
    n.SetZero();

    for(int i = regionBeginIndex; i < (regionBeginIndex+regionSize-1); i++)
        for(int j = i+1; j < (regionBeginIndex+regionSize); j++)
            n += (star[i%star.size()].F()->N() - star[j%star.size()].F()->N());

    return n.SquaredNorm() / std::pow(regionSize, 2);
}


double localCombinatorialEnergy(MyVertex* v,
                                PartitionedVertexStar* partitioning)
{
    std::vector<vcg::face::Pos<MyFace>> star;
    vcg::face::VFOrderedStarFF(vcg::face::Pos<MyFace>(v->VFp(), v), star);
        
    double energy = -1.0;

    if(star.size() <= 3 || v->IsB())
        energy = 0.0;
    else
    {
        // consider all possible cardinalities of a region resulting from a partitioning of the star
        for(int rSize = 2; rSize <= (star.size()-2); rSize++)
        {
            // consider all possible regions having rSize faces and not crossing the vector boundaries.
            // rBegin is the index of the first face in the region
            for(int rBegin = 0; rBegin < (star.size()-rSize); rBegin++)
            {
                double currRegionEnergy =  regionNormalDeviation(star, rBegin,       rSize);
                double otherRegionEnergy = regionNormalDeviation(star, rBegin+rSize, star.size()-rSize);
                
                double currPartitioningEnergy = currRegionEnergy + otherRegionEnergy;

                if(energy < 0 || currPartitioningEnergy < energy)
                {
                    energy = currPartitioningEnergy;
                    if(partitioning != nullptr)
                    {
                        partitioning->star = star;
                        partitioning->rBegin = rBegin;
                        partitioning->rSize = rSize;
                    }
                }
            }
        }
    }

    return energy;
}