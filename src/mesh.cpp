#include "mesh.hpp"


void updateFaceStars(MyMesh& m, StarVertAttrHandle& stars)
{
    int vIndex;
    MyMesh::VertexPointer v;
    std::vector<vcg::face::Pos<MyFace>> currStarPos;

    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
        stars[vIter].clear();

    vcg::tri::UpdateFlags<MyMesh>::VertexClearV(m);
    for(MyMesh::FaceIterator fIter = m.face.begin(); fIter != m.face.end(); fIter++)
    {        
        for(vIndex = 0; vIndex < 3; vIndex++)
        {
            v = fIter->V(vIndex);
            if(v->IsV())
                continue;
            v->SetV();

            vcg::face::VFOrderedStarFF(vcg::face::Pos<MyFace>(&(*fIter), v), currStarPos);
            for(vcg::face::Pos<MyFace> p : currStarPos)
                stars[v].push_back(p.F());
        }
    }
}


void updateNormalsAndAreas(MyMesh& m, AreaFaceAttrHandle& areas)
{
    vcg::tri::UpdateNormal<MyMesh>::PerFace(m);
    double currNorm;
    
    for(MyMesh::FaceIterator fIter = m.face.begin(); fIter != m.face.end(); fIter++)
    {
        currNorm = fIter->N().Norm();
        areas[fIter] = currNorm / 2.0;
        fIter->N().Normalize();
    }
}