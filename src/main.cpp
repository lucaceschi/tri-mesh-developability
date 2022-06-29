#include <iostream>

#include <vcg/complex/algorithms/mesh_to_matrix.h>

#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/export_ply.h>
template<class M> using ImporterOBJ = vcg::tri::io::ImporterOBJ<M>;
template<class M> using ExporterPLY = vcg::tri::io::ExporterPLY<M>;
using TriMask = vcg::tri::io::Mask;

#include "mesh.hpp"


struct PartitionedVertexStar
{
    std::vector<size_t> starI;   // star of face indices
    int rBegin;                  // index of the starting face within starI for the first region
    int rSize;                   // cardinality of the first region
};


double regionNormalDeviation(const Eigen::MatrixXd& normals,
                             std::vector<size_t> starI,
                             int regionBeginIndex,
                             int regionSize)
{   
    Eigen::RowVector3d n = Eigen::Array3d::Zero();

    for(int i = regionBeginIndex; i < (regionBeginIndex + regionSize - 1); i++)
        for(int j = i+1; j < (regionBeginIndex + regionSize); j++)
            n += (normals.row(starI[i%starI.size()]) - normals.row(starI[j%starI.size()]));

    return n.squaredNorm() / std::pow(regionSize, 2);
}


double localCombinatorialEnergy(MyMesh& mesh,
                                size_t v,
                                const Eigen::MatrixX3d& N,
                                PartitionedVertexStar* partitioning = nullptr)
{ 
    MyMesh::VertexPointer vPointer = &mesh.vert[v];
    std::vector<vcg::face::Pos<MyFace>> star;
    std::vector<size_t> starI;
    vcg::face::VFOrderedStarFF(vcg::face::Pos<MyFace>(vPointer->VFp(), vPointer), star);
    for(vcg::face::Pos<MyFace> p : star)
        starI.push_back(vcg::tri::Index(mesh, p.F()));
        
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
                double currRegionEnergy =  regionNormalDeviation(N, starI, rBegin,       rSize);
                double otherRegionEnergy = regionNormalDeviation(N, starI, rBegin+rSize, star.size()-rSize);
                
                double currPartitioningEnergy = currRegionEnergy + otherRegionEnergy;

                if(energy < 0 || currPartitioningEnergy < energy)
                {
                    energy = currPartitioningEnergy;
                    if(partitioning != nullptr)
                    {
                        partitioning->starI = starI;
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
            partitioning->starI = starI;
            partitioning->rBegin = 0;
            partitioning->rSize = 1;
        }
    }
    else return -1.0;

    return energy;
}


void computeNormals(const Eigen::MatrixXd& vert,
                    const Eigen::MatrixXi& faces,
                    Eigen::MatrixXd& normals,
                    Eigen::ArrayXd& areas)
{
    for(size_t fIndex = 0; fIndex < normals.rows(); fIndex++)
    {
        Eigen::Vector3d edgeAvec = vert.row(faces(fIndex, 1)) - vert.row(faces(fIndex, 0));
        Eigen::Vector3d edgeBvec = vert.row(faces(fIndex, 2)) - vert.row(faces(fIndex, 0));

        normals.row(fIndex) = edgeAvec.cross(edgeBvec);
        areas(fIndex) = normals.row(fIndex).norm() / 2.0;
        normals.row(fIndex).normalize();
    }
}


int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        std::cout << "Missing args" << std::endl;
        return 1;
    }
    
    MyMesh m;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd N;
    Eigen::ArrayXd  A;

    int loadMask;
    if(ImporterOBJ<MyMesh>::Open(m, argv[1], loadMask) != ImporterOBJ<MyMesh>::E_NOERROR)
    {
        std::cout << "Error reading input file" << std::endl;
        return 1;
    }
    std::cout << "Loaded file with mask " << loadMask << std::endl;

    vcg::tri::RequireVFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
    vcg::tri::RequireFFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
        
    vcg::tri::MeshToMatrix<MyMesh>::GetTriMeshData(m, F, V);
    N = Eigen::MatrixXd(m.FN(), 3);
    A = Eigen::ArrayXd(m.FN());

    computeNormals(V, F, N, A);

    double totEnergy = 0.0;
    for(size_t v = 0; v < V.rows(); v++)
        totEnergy += localCombinatorialEnergy(m, v, N);
    std::cout << "Total energy = " << totEnergy << std::endl;

    return 0;
}