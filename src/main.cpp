#include <iostream>
#include <iomanip>

#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/export_ply.h>
template<class M> using ImporterOBJ = vcg::tri::io::ImporterOBJ<M>;
template<class M> using ExporterPLY = vcg::tri::io::ExporterPLY<M>;
using TriMask = vcg::tri::io::Mask;

#include "mesh.hpp"
#include "energy.hpp"
#include "energy_grad.hpp"
#include "opt.hpp"
#include "remeshing.hpp"
#define POSTPROCESSING_ANGLE_THRESHOLD 0.3

#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>

#include <time.h>


// -------------------------------------------------------------------------------------------------
// MAIN

int main(int argc, char* argv[])
{
    // >> Mesh loading <<

    MyMesh m;
    double boundingBoxDiag;

    int loadMask;
    if(ImporterOBJ<MyMesh>::Open(m, argv[1], loadMask) != ImporterOBJ<MyMesh>::E_NOERROR)
    {
        std::cout << "Error reading input file" << std::endl;
        return 1;
    }
    std::cout << "Loaded " << argv[1] << " with mask " << loadMask << std::endl;
    
    // >> Mesh preprocessing <<

    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceAdj(m);

    MeshPostProcessing<MyMesh> postProcessing(true, true, POSTPROCESSING_ANGLE_THRESHOLD);
    if(postProcessing.process(m))
        std::cout << "Applied initial remeshing" << std::endl;

    vcg::tri::UpdateBounding<MyMesh>::Box(m);
    boundingBoxDiag = m.bbox.Diag();
    vcg::tri::UpdatePosition<MyMesh>::Scale(m, 1.0 / boundingBoxDiag);

    // >> Optimization <<

    Optimizer* opt;

    // profiling
    double dt;
    std::vector<double> dts;
    clock_t clockStart;

    if(argv[2][0] == 'f')
    {
        int maxFunEval = atoi(argv[3]);
        double eps = atof(argv[4]);
        double stepSize = atof(argv[5]);

        opt = new FixedStepOpt(m, maxFunEval, eps, stepSize);
    }
    else if(argv[2][0] == 'b')
    {
        int maxFunEval = atoi(argv[3]);
        double eps = atof(argv[4]);
        
        opt = new BacktrackingOpt(m, maxFunEval, eps);
    }
    /*else
    {
        int maxFunEval = atoi(argv[3]);
        double eps = atof(argv[4]);
        
        opt = new LewisOvertonOpt(m, maxFunEval, eps);
    }*/

    clockStart = clock();
    while(opt->step())
    {
        opt->printStats();
                
        if(postProcessing.process(m))
        {
            std::cout << "Mesh topology has been altered" << std::endl;
            opt->reset();
        }
        
        dt = ((double)(clock() - clockStart) / CLOCKS_PER_SEC);
        dts.push_back(dt);
        
        clockStart = clock();
    }

    delete opt;

    // compute mean dt
    dt = 0;
    for(double currDt : dts)
        dt += currDt;
    dt /= dts.size();
    std::cout << "Mean dt: " << dt << std::endl;
    
    vcg::tri::UpdatePosition<MyMesh>::Scale(m, boundingBoxDiag);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", TriMask::IOM_NONE, false);

    return 0;
}
