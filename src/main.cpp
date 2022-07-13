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
#include "mesh_matrix.hpp"
#include "energy.hpp"
#include "energy_grad.hpp"

#include <time.h>

#define MIN_ANGLE 0.4


// -------------------------------------------------------------------------------------------------
// MAIN

int main(int argc, char* argv[])
{
    // >> Mesh loading <<

    MyMesh m;

    int loadMask;
    if(ImporterOBJ<MyMesh>::Open(m, argv[1], loadMask) != ImporterOBJ<MyMesh>::E_NOERROR)
    {
        std::cout << "Error reading input file" << std::endl;
        return 1;
    }
    std::cout << "Loaded " << argv[1] << " with mask " << loadMask << std::endl;
    
    // >> Mesh preprocessing <<

    vcg::tri::RequireVFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
    vcg::tri::RequireFFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    vcg::tri::RequirePerVertexFlags<MyMesh>(m);
    vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceAdj(m);

    // >> Mesh 2 matrices convertion <<

    Matrix3Xd V;
    Matrix3Xi F;
    Matrix3Xd N;
    ArrayXd A;
    MatrixXi S;
    ArrayXb B;

    getMeshVF(m, V, F);
    getMeshStars(m, S);
    getMeshBorders(m, B);

    // >> Gradient descent with mesh post processing<<

    Matrix3Xd G(V.rows(), 3);

    int nSteps = atoi(argv[2]);
    float stepSize = atof(argv[3]);
    double totEnergy;

    double dt;
    std::vector<double> dts;
    clock_t clockStart;

    MyMesh::FaceIterator fIter;
    double Ai, Aj, Ak;
    int edgeToFlip;
    bool edgeFlipped;

    for(int step = 0; step < nSteps; step++)
    {
        clockStart = clock();
        computeNormals(V, F, N, A);
        totEnergy = 0.0;

        combinatorialEnergyGrad(V, F, N, A, S, B, G, [&totEnergy](double localEnergy, size_t v) {
            totEnergy += localEnergy;
        });

        V -= (G * stepSize);

        for(size_t v = 0; v < V.rows(); v++)
            m.vert[v].P() = vcg::Point3d(V(v, 0), V(v, 1), V(v, 2));
        
        edgeFlipped = false;
        for(fIter = m.face.begin(); fIter != m.face.end(); fIter++)
        {
            Ai = vcg::face::WedgeAngleRad<MyFace>(*fIter, 0);
            Aj = vcg::face::WedgeAngleRad<MyFace>(*fIter, 1);
            Ak = vcg::face::WedgeAngleRad<MyFace>(*fIter, 2);

            edgeToFlip = -1;
            if(Ai + Aj < MIN_ANGLE)
                edgeToFlip = 0;
            else if(Aj + Ak < MIN_ANGLE)
                edgeToFlip = 1;
            else if(Ak + Ai < MIN_ANGLE)
                edgeToFlip = 2;

            if(edgeToFlip >= 0 && vcg::face::CheckFlipEdge(*fIter, edgeToFlip))
            {
                vcg::face::FlipEdge<MyFace>(*fIter, edgeToFlip);
                edgeFlipped = true;
            }
        }

        if(edgeFlipped)
        {
            std::cout << "Mesh topology has been altered" << std::endl;
            vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
            vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
            vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceAdj(m);
            getMeshVF(m, V, F);
            getMeshStars(m, S);
            getMeshBorders(m, B);
        }

        dt = ((double)(clock() - clockStart) / CLOCKS_PER_SEC);
        std::cout << "[MAT] step #" << step << ": Energy=" << totEnergy << "\tTime=" << dt << std::endl;
        dts.push_back(dt);
    }

    dt = 0;
    for(double currDt : dts)
        dt += currDt;
    dt /= dts.size();
    std::cout << dt << std::endl;

    for(size_t v = 0; v < V.rows(); v++)
        m.vert[v].P() = vcg::Point3d(V(v, 0), V(v, 1), V(v, 2));

    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", TriMask::IOM_NONE, false);

    return 0;
}