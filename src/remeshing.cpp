#include "remeshing.hpp"

#include <vcg/simplex/face/topology.h>
#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/flag.h>

#include "mesh.hpp"


template <class MeshType>
MeshPostProcessing<MeshType>::MeshPostProcessing(bool doEdgeFlipping,
                                                 bool doEdgeCollapsing,
                                                 double angleThreshold) :
    doEdgeFlipping(doEdgeFlipping),
    doEdgeCollapsing(doEdgeCollapsing),
    angleThreshold(angleThreshold)
{}

template MeshPostProcessing<MyMesh>::MeshPostProcessing(bool doEdgeFlipping,
                                                        bool doEdgeCollapsing,
                                                        double angleThreshold);


template <class MeshType>
bool MeshPostProcessing<MeshType>::process(MeshType& m)
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FaceIterator FaceIterator;
    
    FaceIterator fIter;
    double faceAngles[3];
    int edgeToFlip;
    int edgeToCollapse;
    bool meshModified = false;
    
    for(fIter = m.face.begin(); fIter != m.face.end(); fIter++)
    {
        if(fIter->IsD())
            continue;
        
        faceAngles[0] = vcg::face::WedgeAngleRad<FaceType>(*fIter, 0);
        faceAngles[1] = vcg::face::WedgeAngleRad<FaceType>(*fIter, 1);
        faceAngles[2] = vcg::face::WedgeAngleRad<FaceType>(*fIter, 2);

        edgeToFlip = -1;
        edgeToCollapse = -1;

        for(int i = 0; i < 3; i++)
            if(faceAngles[i] < angleThreshold)
            {
                if(faceAngles[(i+1)%3] < angleThreshold)
                {
                    edgeToFlip = i;
                    break;
                }
                else if(faceAngles[(i+2)%3] >= angleThreshold)
                {
                    edgeToCollapse = (i+1)%3;
                    break;
                }
            }

        if(doEdgeFlipping && edgeToFlip >= 0 && vcg::face::CheckFlipEdge(*fIter, edgeToFlip))
        {
            vcg::face::FlipEdge<FaceType>(*fIter, edgeToFlip);
            meshModified = true;
        }
        else if(doEdgeCollapsing && edgeToCollapse >= 0 && vcg::face::FFLinkCondition(*fIter, edgeToCollapse))
        {
            vcg::face::FFEdgeCollapse<MeshType>(m, *fIter, edgeToCollapse);
            meshModified = true;
        }
    }

    if(meshModified)
    {
        vcg::tri::Allocator<MeshType>::CompactFaceVector(m);
        vcg::tri::Allocator<MeshType>::CompactVertexVector(m);
        vcg::tri::UpdateTopology<MeshType>::VertexFace(m);
        vcg::tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(m);
    }

    return meshModified;
}

template bool MeshPostProcessing<MyMesh>::process(MyMesh& m);