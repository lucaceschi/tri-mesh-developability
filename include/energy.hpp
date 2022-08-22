#ifndef ENERGY_HPP
#define ENERGY_HPP

#include "mesh.hpp"


struct StarPartitioning
{
    Star* star;
    int rBegin;   // index of the starting face within star for the first region
    int rSize;    // cardinality of the first region
};

double combinatorialEnergy(MyMesh& m,
                           StarVertAttrHandle& vAttrStar);

double localCombinatorialEnergy(MyMesh::VertexPointer v,
                                MyMesh& m,
                                StarVertAttrHandle& vAttrStar,
                                StarPartitioning* partitioning = nullptr);

double regionNormalDeviation(const StarPartitioning& partitioning,
                             bool region,
                             MyMesh& m);


#endif