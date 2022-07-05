#ifndef ENERGY_HPP
#define ENERGY_HPP

#include <vector>

#include "mesh.hpp"


struct PartitionedVertexStar
{
    std::vector<vcg::face::Pos<MyFace>> star;   // star of faces
    int rBegin;                                 // index of the first face of the first region
    int rSize;                                  // cardinality of the first region
};


double regionNormalDeviation(std::vector<vcg::face::Pos<MyFace>> star,
                             int regionBeginIndex,
                             int regionSize);

double localCombinatorialEnergy(MyVertex* v,
                                PartitionedVertexStar* partitioning = nullptr);

#endif