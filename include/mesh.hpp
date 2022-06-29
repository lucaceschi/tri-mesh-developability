#ifndef MESH_HPP
#define MESH_HPP

#include <vcg/complex/complex.h>

class MyVertex;
class MyEdge;
class MyFace;

struct MyUsedTypes : public vcg::UsedTypes< vcg::Use<MyVertex>::AsVertexType,
                                            vcg::Use<MyEdge>::AsEdgeType,
                                            vcg::Use<MyFace>::AsFaceType >{};
                                           
class MyVertex : public vcg::Vertex< MyUsedTypes,
                                     vcg::vertex::Coord3d,
                                     vcg::vertex::VFAdj,
                                     vcg::vertex::BitFlags >{};

class MyFace : public vcg::Face< MyUsedTypes,
                                 vcg::face::VertexRef,
                                 vcg::face::VFAdj,
                                 vcg::face::FFAdj,
                                 vcg::face::Normal3d,
                                 vcg::face::BitFlags >{};

class MyEdge : public vcg::Edge< MyUsedTypes >{};

class MyMesh : public vcg::tri::TriMesh< std::vector<MyVertex>,
                                         std::vector<MyFace>,
                                         std::vector<MyEdge> >{};

#endif