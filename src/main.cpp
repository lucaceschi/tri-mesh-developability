#include <iostream>
#include <iomanip>

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
#include "opt.hpp"
#include "remeshing.hpp"

#include <time.h>

#define POSTPROCESSING_ANGLE_THRESHOLD 0.3


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

    Matrix3Xd V(m.VN(), 3);
    Matrix3Xi F(m.FN(), 3);
    Matrix3Xd N(m.FN(), 3);
    ArrayXd A(m.FN());
    MatrixXi S;
    ArrayXb B(m.VN());

    getMeshVF(m, V, F);
    getMeshStars(m, S);
    getMeshBorders(m, B);

    // >> Optimization <<

    Optimizer* opt;
    MeshPostProcessing<MyMesh> postProcessing(true, true, POSTPROCESSING_ANGLE_THRESHOLD);

    // profiling
    double dt;
    std::vector<double> dts;
    clock_t clockStart;

    if(argv[2][0] == 'f')
    {
        int maxFunEval = atoi(argv[3]);
        double eps = atof(argv[4]);
        double stepSize = atof(argv[5]);

        opt = new FixedStepOpt(V.rows(), maxFunEval, eps, stepSize);
    }
    else if(argv[2][0] == 'b')
    {
        int maxFunEval = atoi(argv[3]);
        double eps = atof(argv[4]);
        
        opt = new BacktrackingOpt(V.rows(), maxFunEval, eps);
    }
    else
    {
        int maxFunEval = atoi(argv[3]);
        double eps = atof(argv[4]);
        
        opt = new LewisOvertonOpt(V.rows(), maxFunEval, eps);
    }

    clockStart = clock();
    while(opt->step(V, F, N, A, S, B))
    {
        for(size_t v = 0; v < V.rows(); v++)
            m.vert[v].P() = vcg::Point3d(V(v, 0), V(v, 1), V(v, 2));
        
        if(postProcessing.process(m))
        {
            std::cout << "Mesh topology has been altered" << std::endl;
            getMeshVF(m, V, F);
            getMeshStars(m, S);
            getMeshBorders(m, B);
        }
        
        dt = ((double)(clock() - clockStart) / CLOCKS_PER_SEC);
        dts.push_back(dt);
        opt->printStats();
        clockStart = clock();
    }

    delete opt;

    // compute mean dt
    dt = 0;
    for(double currDt : dts)
        dt += currDt;
    dt /= dts.size();
    std::cout << "Mean dt: " << dt << std::endl;

    for(size_t v = 0; v < V.rows(); v++)
        m.vert[v].P() = vcg::Point3d(V(v, 0), V(v, 1), V(v, 2));

    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", TriMask::IOM_NONE, false);

    return 0;
}