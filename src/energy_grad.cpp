#include "energy_grad.hpp"


vcg::Matrix33d faceNormalGrad(MyMesh::FacePointer f,
                               int vIndex,
                               AreaFaceAttrHandle& fAttrArea)
{    
    vcg::Point3d oppositeEdge = f->V2(vIndex)->P() - f->V1(vIndex)->P();
    vcg::Matrix33d grad;

    grad.ExternalProduct(oppositeEdge ^ f->N(), f->N());
    grad /= fAttrArea[f];

    return grad;
}


void regionNormalDeviationGrad(MyMesh::VertexPointer v,
                               StarPartitioning& partitioning,
                               bool region,
                               MyMesh& m,
                               AreaFaceAttrHandle& fAttrArea,
                               StarVertAttrHandle& vAttrStar,
                               GradientVertAttrHandle& vAttrGrad)
{  
    int rBegin = region ? (partitioning.rBegin + partitioning.rSize)       : partitioning.rBegin;
    int rSize  = region ? (partitioning.star->size() - partitioning.rSize) : partitioning.rSize;
    int starSize = partitioning.star->size();
    
    MyMesh::FacePointer faceA, faceB;
    vcg::Point3d normalDiff;

    int faceA_vIndex, faceB_vIndex;
    MyMesh::VertexPointer vertA, vertB;
    
    for(int i = rBegin; i < (rBegin + rSize - 1); i++)
        for(int j = i+1; j < (rBegin + rSize); j++)
        {
            faceA = vAttrStar[v][i % starSize];
            faceB = vAttrStar[v][j % starSize];

            normalDiff = faceA->N() - faceB->N();

            for(faceA_vIndex = 0; faceA_vIndex < 3; faceA_vIndex++)
            {
                vertA = faceA->V(faceA_vIndex);
                vAttrGrad[vertA] += (faceNormalGrad(faceA, faceA_vIndex, fAttrArea).transpose() * normalDiff * 2/std::pow(rSize, 2));
            }

            for(faceB_vIndex = 0; faceB_vIndex < 3; faceB_vIndex++)
            {
                vertB = faceB->V(faceB_vIndex);
                vAttrGrad[vertB] -= (faceNormalGrad(faceB, faceB_vIndex, fAttrArea).transpose() * normalDiff* 2/std::pow(rSize, 2));
            }            
        }
}


double combinatorialEnergyGrad(MyMesh& m,
                               AreaFaceAttrHandle& fAttrArea,
                               StarVertAttrHandle& vAttrStar,
                               GradientVertAttrHandle& vAttrGrad)
{
    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
        vAttrGrad[vIter].SetZero();

    double currEnergy;
    double totEnergy = 0.0;
    StarPartitioning currPart;
    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
    {
        currEnergy = localCombinatorialEnergy(&(*vIter), m, vAttrStar, &currPart);
        totEnergy += currEnergy;

        if(currPart.star->size() <= 3 || vIter->IsB())
            continue;

        regionNormalDeviationGrad(&(*vIter), currPart, 0, m, fAttrArea, vAttrStar, vAttrGrad);
        regionNormalDeviationGrad(&(*vIter), currPart, 1, m, fAttrArea, vAttrStar, vAttrGrad);
    }

    return totEnergy;
}