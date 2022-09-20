#ifndef ENERGY_HPP
#define ENERGY_HPP

#include "mesh.hpp"


struct StarPartitioning
{
    Star* star;
    int rBegin;   // index of the starting face within star for the first region
    int rSize;    // cardinality of the first region
};


/*
 * Compute the combinatorial energy of a mesh
 * as defined in section 4.1.1 of the paper "Developability of Triangle Meshes" by Stein et al. 2018
 */
double combinatorialEnergy(MyMesh& m,
                           StarVertAttrHandle& vAttrStar);


/*
 * Compute the combinatorial energy of a vertex star
 * as defined in section 4.1.1 of the paper "Developability of Triangle Meshes" by Stein et al. 2018
 */
double localCombinatorialEnergy(MyMesh::VertexPointer v,
                                MyMesh& m,
                                StarVertAttrHandle& vAttrStar,
                                StarPartitioning* partitioning = nullptr);


/*
 * Compute the normal deviation of a region of a vertex star
 * as defined in equation 1 of the paper "Developability of Triangle Meshes" by Stein et al. 2018
 */
double regionNormalDeviation(const StarPartitioning& partitioning,
                             bool region,
                             MyMesh& m);


#endif