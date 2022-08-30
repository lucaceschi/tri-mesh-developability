#ifndef REMESHING_HPP
#define REMESHING_HPP


/*
 * To improve numerical stability of the energy optimization, meshes should not present
 * small interior angles;
 * this class implements a post-processing phase that can be applied between iterations
 * in order to get rid of these angles by means of edge flipping and/or edge collapsing
 */
template <class MeshType>
class MeshPostProcessing
{
public:
    MeshPostProcessing(bool doEdgeFlipping,
                       bool doEdgeCollapsing,
                       double angleThreshold);

    /*
     * Perform the post-processing on a mesh
     */
    bool process(MeshType& m);

private:
    bool doEdgeFlipping;
    bool doEdgeCollapsing;
    double angleThreshold;
};


#endif