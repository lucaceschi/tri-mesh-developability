#include <iostream>

#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/space/color4.h>
#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/export_ply.h>

#include "mesh.hpp"


double regionNormalDeviation(std::vector<vcg::face::Pos<MyFace>> star, int regionBeginIndex, int regionSize)
{   
    vcg::Point3d n;
    n.SetZero();

    for(int i = regionBeginIndex; i < (regionBeginIndex+regionSize-1); i++)
        for(int j = i+1; j < (regionBeginIndex+regionSize); j++)
            n += (star[i%star.size()].F()->N() - star[j%star.size()].F()->N()); // TODO migliorare

    return n.SquaredNorm() / std::pow(regionSize, 2);
}


double localCombinatorialEnergy(MyVertex* v)
{
    std::vector<vcg::face::Pos<MyFace>> star;
    vcg::face::VFOrderedStarFF(vcg::face::Pos<MyFace>(v->VFp(), v), star);
        
    // TODO: ritoccare gestione casi speciali |VF|=3 e |VF|<2
    double energy = -1.0;

    if(star.size() == 2)
    {
        energy = 0.0;
    }
    else if(star.size() > 3)
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
                    energy = currPartitioningEnergy;
            }
        }
    }

    return energy;
}


int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::cout << "No input file provided" << std::endl;
        return 1;
    }
    
    MyMesh m;

    //vcg::tri::Hexahedron(m);

    int loadMask;
    if(vcg::tri::io::ImporterOBJ<MyMesh>::Open(m, argv[1], loadMask) != vcg::tri::io::ImporterOBJ<MyMesh>::E_NOERROR)
    {
        std::cout << "Error reading input file" << std::endl;
        return 1;
    }
    std::cout << "Loaded file with mask " << loadMask << std::endl;

    vcg::tri::RequireVFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
    vcg::tri::RequireFFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    
    vcg::tri::RequirePerFaceNormal(m);
    vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);

    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
        (*vIter).Q() = localCombinatorialEnergy(&*vIter);

    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", vcg::tri::io::Mask::IOM_VERTQUALITY, false);

    return 0;
}