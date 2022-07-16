#ifndef REMESHING_HPP
#define REMESHING_HPP


template <class MeshType>
class MeshPostProcessing
{
public:
    MeshPostProcessing(bool doEdgeFlipping,
                       bool doEdgeCollapsing,
                       double angleThreshold);

    bool process(MeshType& m);

private:
    bool doEdgeFlipping;
    bool doEdgeCollapsing;
    double angleThreshold;
};


#endif