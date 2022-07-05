#include "energy_grad.hpp"


vcg::Matrix33d faceNormalGrad(MyFace* f,
                              int v,
                              MyMesh::PerFaceAttributeHandle<double>& triArea)
{
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
}


void regionNormalDeviationGrad(MyVertex* v,
                               std::vector<vcg::face::Pos<MyFace>> star,
                               int regionBeginIndex,
                               int regionSize,
                               MyMesh::PerFaceAttributeHandle<double>& triArea,
                               MyMesh::PerVertexAttributeHandle<vcg::Point3d>& energyGrad)
{
    MyMesh::FacePointer faceA, faceB;
    vcg::Point3d faceNormDiff;
    int faceAvIndex, faceBvIndex;
    MyVertex *faceAv, *faceBv;
    bool sharedVertex;
    
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
                        energyGrad[faceAv] += ((faceNormalGrad(faceA, faceAvIndex, triArea) - faceNormalGrad(faceB, faceBvIndex, triArea)).transpose() * faceNormDiff * 2/std::pow(regionSize, 2));
                        sharedVertex = true;
                        break;
                    }
                }

                if(!sharedVertex)
                    energyGrad[faceAv] += (faceNormalGrad(faceA, faceAvIndex, triArea).transpose() * faceNormDiff * 2/std::pow(regionSize, 2));
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
                    energyGrad[faceBv] -= (faceNormalGrad(faceB, faceBvIndex, triArea).transpose() * faceNormDiff * 2/std::pow(regionSize, 2));
            }
        }
}


void combinatorialEnergyGrad(MyMesh& mesh,
                             MyMesh::PerVertexAttributeHandle<vcg::Point3d>& energyGrad,
                             std::function<void(double, MyVertex*)> localEnergyCallback)
{
    MyMesh::PerFaceAttributeHandle<double> triArea = vcg::tri::Allocator<MyMesh>::GetPerFaceAttribute<double>(mesh, "Area");

    PartitionedVertexStar currPart;
    double currEnergy;

    // set to zero the area of each face
    for(MyMesh::FaceIterator fIter = mesh.face.begin(); fIter != mesh.face.end(); fIter++)
        triArea[fIter] = 0.0;

    for(MyMesh::VertexIterator vIter = mesh.vert.begin(); vIter != mesh.vert.end(); vIter++)
    {
        currEnergy = localCombinatorialEnergy(&*vIter, &currPart);
        localEnergyCallback(currEnergy, &*vIter);

        if(currPart.star.size() <= 3 || vIter->IsB())
            continue;

        regionNormalDeviationGrad(
            &*vIter, currPart.star, currPart.rBegin, currPart.rSize,
            triArea, energyGrad
        );

        regionNormalDeviationGrad(
            &*vIter, currPart.star, currPart.rBegin+currPart.rSize, currPart.star.size() - currPart.rSize,
            triArea, energyGrad
        );
    }
}