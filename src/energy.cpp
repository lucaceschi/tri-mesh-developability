#include "energy.hpp"


double combinatorialEnergy(MyMesh& m,
                           StarVertAttrHandle& vAttrStar)
{
    double totEnergy = 0.0;

    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
        totEnergy += localCombinatorialEnergy(&(*vIter), m, vAttrStar);

    return totEnergy;
}


double localCombinatorialEnergy(MyMesh::VertexPointer v,
                                MyMesh& m,
                                StarVertAttrHandle& vAttrStar,
                                StarPartitioning* outP)
{
    double energy = -1.0;
    StarPartitioning currP = {
        &(vAttrStar[v]),    // star of faces around v
        0,                  // rBegin
        0                   // rSize
    };

    if(outP != nullptr)
        outP->star = currP.star;

    if(currP.star->size() <= 3 || v->IsB())
        energy = 0.0;
    else
    {
        // consider all possible cardinalities of a region resulting from a partitioning of the star
        for(currP.rSize = 2; currP.rSize <= (currP.star->size() - 2); currP.rSize++)
        {
            // consider all possible regions having rSize faces and not crossing the vector boundaries.
            // rBegin is the index of the first face in the region
            for(currP.rBegin = 0; currP.rBegin < (currP.star->size() - currP.rSize); currP.rBegin++)
            {                
                double currRegionEnergy  = regionNormalDeviation(currP, 0, m);
                double otherRegionEnergy = regionNormalDeviation(currP, 1, m);                
                double currPartitioningEnergy = currRegionEnergy + otherRegionEnergy;

                if(energy < 0 || currPartitioningEnergy < energy)
                {
                    energy = currPartitioningEnergy;
                    if(outP != nullptr)
                    {
                        outP->rBegin = currP.rBegin;
                        outP->rSize = currP.rSize;
                    }
                }
            }
        }
    }

    return energy;
}


double regionNormalDeviation(const StarPartitioning& partitioning,
                             bool region,
                             MyMesh& m)
{
    int rBegin = region ? (partitioning.rBegin + partitioning.rSize)       : partitioning.rBegin;
    int rSize  = region ? (partitioning.star->size() - partitioning.rSize) : partitioning.rSize;
    int starSize = partitioning.star->size();

    vcg::Point3d normDiff;
    normDiff.SetZero();

    double regionNormalDev = 0.0;

    for(int i = rBegin; i < (rBegin + rSize - 1); i++)
        for(int j = i+1; j < (rBegin + rSize); j++)
        {
            normDiff = ( partitioning.star->at(i % starSize)->N() - partitioning.star->at(j % starSize)->N() );
            regionNormalDev += ( normDiff.SquaredNorm() / std::pow(rSize, 2) );
        }

    return regionNormalDev;
}
