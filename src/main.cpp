#include <iostream>

#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/space/color4.h>
#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_ply.h>
template<class M> using ImporterOBJ = vcg::tri::io::ImporterOBJ<M>;
template<class M> using ExporterPLY = vcg::tri::io::ExporterPLY<M>;
using TriMask = vcg::tri::io::Mask;

#include "mesh.hpp"


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


struct PartitionedVertexStar
{
    std::vector<vcg::face::Pos<MyFace>> star;   // star of faces
    int rBegin;                                 // index of the first face of the first region
    int rSize;                                  // cardinality of the first region
};


double localCombinatorialEnergy(MyVertex* v,
                                PartitionedVertexStar* partitioning = nullptr)
{
    std::vector<vcg::face::Pos<MyFace>> star;
    vcg::face::VFOrderedStarFF(vcg::face::Pos<MyFace>(v->VFp(), v), star);
        
    double energy = -1.0;

    if(star.size() > 3)
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
    else if(star.size() == 2)
    {
        energy = 0.0;
        if(partitioning != nullptr)
        {
            partitioning->star = star;
            partitioning->rBegin = 0;
            partitioning->rSize = 1;
        }
    }
    else return -1.0;

    return energy;
}


void combinatorialEnergyGrad(MyMesh& mesh,
                             MyMesh::PerVertexAttributeHandle<vcg::Point3d> &energyGrad,
                             std::function<void(double, MyVertex*)> localEnergyCallback)
{
    MyMesh::PerFaceAttributeHandle<double> triArea = vcg::tri::Allocator<MyMesh>::GetPerFaceAttribute<double>(mesh, "Area");

    PartitionedVertexStar currPart;
    double currEnergy;

    MyMesh::FacePointer faceA, faceB;
    vcg::Point3d faceNormDiff;
    int vIndex;

    // set to zero the area of each face
    for(MyMesh::FaceIterator fIter = mesh.face.begin(); fIter != mesh.face.end(); fIter++)
        triArea[fIter] = 0.0;

    auto faceNormalGrad = [&triArea](MyFace* f, int v) -> vcg::Matrix33d {
        vcg::Point3d oppositeEdge = f->V1(v)->P() - f->V2(v)->P();
        vcg::Matrix33d grad;

        // if unknown, compute&store triangle area
        if(triArea[f] == 0.0)
        {
            vcg::Point3d properEdge = f->V(v)->P() - f->V1(v)->P();
            triArea[f] = (oppositeEdge ^ properEdge).Norm() / 2.0;
        }

        grad.ExternalProduct(oppositeEdge ^ f->N(), f->N());
        grad /= triArea[f];

        return grad;
    };

    auto regionNormalDeviationGrad = [&](MyVertex* v, std::vector<vcg::face::Pos<MyFace>> star, int regionBeginIndex, int regionSize) {
        for(int i = regionBeginIndex; i < (regionBeginIndex+regionSize-1); i++)
            for(int j = i+1; j < (regionBeginIndex+regionSize); j++)
            {
                faceA = star[i%star.size()].F();
                faceB = star[j%star.size()].F();

                faceNormDiff = faceA->N() - faceB->N();

                for(vIndex = 0; vIndex < 3; vIndex++)
                {
                    if(faceA->V(vIndex) == v)
                        energyGrad[v] += ((faceNormalGrad(faceA, vIndex) - faceNormalGrad(faceB, vIndex)).transpose() * faceNormDiff * 2);
                    else
                        energyGrad[v] += (faceNormalGrad(faceA, vIndex).transpose() * faceNormDiff * 2);
                }

                for(vIndex = 0; vIndex < 3; vIndex++)
                {
                    if(faceB->V(vIndex) != v)
                        energyGrad[v] += (faceNormalGrad(faceB, vIndex).transpose() * faceNormDiff * 2);
                }
            }
    };

    for(MyMesh::VertexIterator vIter = mesh.vert.begin(); vIter != mesh.vert.end(); vIter++)
    {
        currEnergy = localCombinatorialEnergy(&*vIter, &currPart);
        localEnergyCallback(currEnergy, &*vIter);

        if(currPart.star.size() == 3)
            continue;

        regionNormalDeviationGrad(&*vIter, currPart.star, currPart.rBegin, currPart.rSize);
        regionNormalDeviationGrad(&*vIter, currPart.star, currPart.rBegin+currPart.rSize, currPart.star.size() - currPart.rSize);

        energyGrad[vIter] /= std::pow(currPart.star.size(), 2);
    }
}


int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        std::cout << "Missing args" << std::endl;
        return 1;
    }
    
    MyMesh m;

    int loadMask;
    if(ImporterOBJ<MyMesh>::Open(m, argv[1], loadMask) != ImporterOBJ<MyMesh>::E_NOERROR)
    {
        std::cout << "Error reading input file" << std::endl;
        return 1;
    }
    std::cout << "Loaded " << argv[1] << " with mask " << loadMask << std::endl;

    vcg::tri::RequireVFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
    vcg::tri::RequireFFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    
    vcg::tri::RequirePerFaceNormal(m);

    MyMesh::PerVertexAttributeHandle<vcg::Point3d> energyGrad = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<vcg::Point3d>(m, "Combinatorial energy gradient");
    double totEnergy;

    // perform standard gradient descent with argv[2] steps and a stepsize of argv[3]
    int nSteps = atoi(argv[2]);
    float stepSize = atof(argv[3]);
    for(int step = 0; step < nSteps; step++)
    {
        vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);
        totEnergy = 0.0;

        // reset gradient of each vertex
        for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
            energyGrad[vIter].SetZero();

        combinatorialEnergyGrad(m, energyGrad, [&totEnergy](double energy, MyVertex* v) {
            totEnergy += energy;
        });

        std::cout << "Step #" << step << "\tEnergy=" << totEnergy << std::endl;

        for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
            (*vIter).P() += (energyGrad[vIter] * stepSize);
    }

    totEnergy = 0.0;
    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
        totEnergy += localCombinatorialEnergy(&*vIter);

    std::cout << "Final Energy= " << totEnergy << std::endl;

    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", TriMask::IOM_NONE, false);

    return 0;
}