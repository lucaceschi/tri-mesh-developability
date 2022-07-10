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

   // >> Gradient <<

    Matrix3Xd G(V.rows(), 3);

    int nSteps = atoi(argv[2]);
    float stepSize = atof(argv[3]);
    double totEnergy;

    double dt;
    std::vector<double> dts;
    clock_t clockStart;

    for(int step = 0; step < nSteps; step++)
    {
        clockStart = clock();
        computeNormals(V, F, N, A);
        totEnergy = 0.0;

        combinatorialEnergyGrad(V, F, N, A, S, B, G, [&totEnergy](double localEnergy, size_t v) {
            totEnergy += localEnergy;
        });

        V -= (G * stepSize);

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