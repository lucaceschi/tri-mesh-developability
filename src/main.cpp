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


void combinatorialEnergyGrad(MyMesh& mesh,
                             MyMesh::PerVertexAttributeHandle<vcg::Point3d> &energyGrad,
                             std::function<void(double, MyVertex*)> localEnergyCallback)
{
    MyMesh::PerFaceAttributeHandle<double> triArea = vcg::tri::Allocator<MyMesh>::GetPerFaceAttribute<double>(mesh, "Area");

    PartitionedVertexStar currPart;
    double currEnergy;

    MyMesh::FacePointer faceA, faceB;
    vcg::Point3d faceNormDiff;
    int faceAvIndex, faceBvIndex;
    MyVertex *faceAv, *faceBv;
    bool sharedVertex;

    // set to zero the area of each face
    for(MyMesh::FaceIterator fIter = mesh.face.begin(); fIter != mesh.face.end(); fIter++)
        triArea[fIter] = 0.0;

    auto faceNormalGrad = [&triArea](MyFace* f, int v) -> vcg::Matrix33d {
        vcg::Point3d oppositeEdge = f->V2(v)->P() - f->V1(v)->P();
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

                for(faceAvIndex = 0; faceAvIndex < 3; faceAvIndex++)
                {
                    faceAv = faceA->V(faceAvIndex);
                    
                    // check if the current vertex of faceA is shared with faceB
                    sharedVertex = false;
                    for(faceBvIndex = 0; faceBvIndex < 3; faceBvIndex++)
                    {
                        faceBv = faceB->V(faceBvIndex);
                        if(faceAv == faceBv)
                        {
                            energyGrad[faceAv] += ((faceNormalGrad(faceA, faceAvIndex) - faceNormalGrad(faceB, faceBvIndex)).transpose() * faceNormDiff * 2/std::pow(regionSize, 2));
                            sharedVertex = true;
                            break;
                        }
                    }

                    if(!sharedVertex)
                        energyGrad[faceAv] += (faceNormalGrad(faceA, faceAvIndex).transpose() * faceNormDiff * 2/std::pow(regionSize, 2));
                }

                for(faceBvIndex = 0; faceBvIndex < 3; faceBvIndex++)
                {
                    faceBv = faceB->V(faceBvIndex);

                    // make sure the current vertex of faceB is not shared with faceA
                    sharedVertex = false;
                    for(faceAvIndex = 0; faceAvIndex < 3; faceAvIndex++)
                    {
                        faceAv = faceA->V(faceAvIndex);
                        if(faceAv == faceBv)
                        {
                            sharedVertex = true;
                            break;
                        }
                    }

                    if(!sharedVertex)
                        energyGrad[faceBv] -= (faceNormalGrad(faceB, faceBvIndex).transpose() * faceNormDiff * 2/std::pow(regionSize, 2));
                }
            }
    };

    for(MyMesh::VertexIterator vIter = mesh.vert.begin(); vIter != mesh.vert.end(); vIter++)
    {
        currEnergy = localCombinatorialEnergy(&*vIter, &currPart);
        localEnergyCallback(currEnergy, &*vIter);

        if(currPart.star.size() <= 3 || vIter->IsB())
            continue;

        regionNormalDeviationGrad(&*vIter, currPart.star, currPart.rBegin, currPart.rSize);
        regionNormalDeviationGrad(&*vIter, currPart.star, currPart.rBegin+currPart.rSize, currPart.star.size() - currPart.rSize);
    }
}


int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        std::cout << "Missing args" << std::endl;
        return 1;
    }
    
    // >> Mesh loading & preprocessing <<

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
    vcg::tri::RequirePerVertexFlags<MyMesh>(m);
    vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceAdj(m);
    
    vcg::tri::RequirePerFaceNormal(m);

    // >> Gradient <<

    MyMesh::PerVertexAttributeHandle<vcg::Point3d> energyGrad = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<vcg::Point3d>(m, "Combinatorial energy gradient");
    double totEnergy;

    double dt;
    std::vector<double> dts;
    clock_t clockStart;

    // perform simple gradient descent with argv[2] steps and a stepsize of argv[3]
    int nSteps = atoi(argv[2]);
    float stepSize = atof(argv[3]);
    for(int step = 0; step < nSteps; step++)
    {
        clockStart = clock();
        vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);
        totEnergy = 0.0;

        // reset gradient of each vertex
        for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
            energyGrad[vIter].SetZero();

        combinatorialEnergyGrad(m, energyGrad, [&totEnergy](double energy, MyVertex* v) {
            totEnergy += energy;
        });

        for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
            (*vIter).P() -= (energyGrad[vIter] * stepSize);

        dt = ((double)(clock() - clockStart) / CLOCKS_PER_SEC);
        std::cout << "Step #" << step << ": Energy=" << totEnergy << "\tTime=" << dt << std::endl;
        dts.push_back(dt);    
    }

    dt = 0;
    for(double currDt : dts)
        dt += currDt;
    dt /= dts.size();
    std::cout << "Mean time per step: " << dt << std::endl;

    totEnergy = 0.0;
    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
        totEnergy += localCombinatorialEnergy(&*vIter);

    std::cout << "Final energy: " << totEnergy << std::endl;

    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", TriMask::IOM_NONE, false);

    return 0;
}