#include <stdlib.h>
#include <iostream>
#include <string>

#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_ply.h>
template<class M> using ImporterOBJ = vcg::tri::io::ImporterOBJ<M>;
template<class M> using ExporterPLY = vcg::tri::io::ExporterPLY<M>;
using TriMask = vcg::tri::io::Mask;

#include "mesh.hpp"
#include "energy.hpp"
#include "energy_grad.hpp"
#include "opt.hpp"
#include "remeshing.hpp"
#define POSTPROCESSING_ANGLE_THRESHOLD_RAD 0.3

#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/position.h>

#include <time.h>


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// MAIN

int main(int argc, char* argv[])
{
    /*/ ARGS PARSING /*/

    if(argc < 6 || (argv[2][0] != 'f' && argv[2][0] != 'b') || (argv[2][0] == 'b' && argc < 7))
    {
        std::cerr << "Wrong argument count!" << std::endl;
        std::cout << "Usage: exec "
                     "[path to mesh] "              // argv[1]
                     "[optimizer] "                 // argv[2]
                     "[max no. fun evals] "         // argv[3]
                     "[stop gradient threshold] "   // argv[4]
                     "[step size] "                 // argv[5]
                     "{tau}"                        // argv[6]
                  << std::endl;
        std::cout << "optimizer can be: " << std::endl
                  << "\t f = Gradient method opt with fixed step size" << std::endl
                  << "\t b = Gradient method opt with backtracking line search (Armijo condition)" << std::endl
                  << std::endl;
        std::cout << "Optimization will stop when the max number of function evaluations is reached, "
                     "or the squared gradient goes below the given threshold." << std::endl
                  << std::endl;
        std::cout << "When using backtracking opt, step size is progressively scaled down "
                  << "by factors of tau according to the Armijo condition" << std::endl
                  << std::endl;
        std::cout << "Example usage:" << std::endl
                  << "\t develop ../assets/bunny.obj b 400 0 0.0001 0.8" << std::endl;

        exit(EXIT_FAILURE);
    }

    int maxFunEval = 0;
    double eps = 0.0;
    double stepSize = 0.0;
    double tau = 0.0;

    try
    {
        maxFunEval = std::stoi(std::string(argv[3]));
        eps        = std::stod(std::string(argv[4]));
        stepSize   = std::stod(std::string(argv[5]));
        if(argv[2][0] == 'b')
            tau = std::stod(std::string(argv[6]));
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    
    /*/ MESH LOADING /*/

    MyMesh m;
    double boundingBoxDiag;

    int loadMask;
    int err;

    err = ImporterOBJ<MyMesh>::Open(m, argv[1], loadMask);
    if(err != ImporterOBJ<MyMesh>::E_NOERROR)
    {
        std::cerr << "Error reading input file (ImporterOBJ::OBJError " << err << ")" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "Loaded " << argv[1] << " with mask " << loadMask << std::endl;
    
    /*/ MESH PREPROCESSING /*/

    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceAdj(m);

    for(MyMesh::FaceIterator fIter = m.face.begin(); fIter != m.face.end(); fIter++)
        for(int i = 0; i < 3; i++)
            if(!vcg::face::IsManifold(*fIter, i))
            {
                std::cerr << "Input mesh must be manifold" << std::endl;
                exit(EXIT_FAILURE);
            }

    // perform an initial remeshing if necessary
    MeshPostProcessing<MyMesh> postProcessing(true, true, POSTPROCESSING_ANGLE_THRESHOLD_RAD);
    if(postProcessing.process(m))
        std::cout << "An initial remeshing has been applied" << std::endl;

    vcg::tri::UpdateBounding<MyMesh>::Box(m);
    boundingBoxDiag = m.bbox.Diag();
    vcg::tri::UpdatePosition<MyMesh>::Scale(m, 1.0 / boundingBoxDiag);

    /*/ OPTIMIZATION /*/

    Optimizer* opt;

    // variables for profiling
    double dt;
    std::vector<double> dts;
    clock_t clockStart;

    if(argv[2][0] == 'f')
        opt = new FixedStepOpt(m, maxFunEval, eps, stepSize);
    else
        opt = new BacktrackingOpt(m, maxFunEval, eps, stepSize, tau);

    clockStart = clock();
    while(opt->step())
    {
        opt->printStats();
                
        if(postProcessing.process(m))
        {
            std::cout << "Remeshing applied" << std::endl;
            opt->reset();
        }
        
        dt = ((double)(clock() - clockStart) / CLOCKS_PER_SEC);
        dts.push_back(dt);
        
        clockStart = clock();
    }

    delete opt;

    // compute mean delta time
    dt = 0;
    for(double currDt : dts)
        dt += currDt;
    dt /= dts.size();
    std::cout << "Mean delta time: " << dt << std::endl;

    vcg::tri::UpdatePosition<MyMesh>::Scale(m, boundingBoxDiag);
    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", TriMask::IOM_NONE, false);

    return 0;
}
