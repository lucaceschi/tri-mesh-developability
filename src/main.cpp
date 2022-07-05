#include <iostream>

#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/space/color4.h>
#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/import_obj.h>
#include <wrap/io_trimesh/export_ply.h>
template<class M> using ImporterOBJ = vcg::tri::io::ImporterOBJ<M>;
template<class M> using ExporterPLY = vcg::tri::io::ExporterPLY<M>;
using TriMask = vcg::tri::io::Mask;

#include "mesh.hpp"
#include "energy_grad.hpp"


int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        std::cout << "Missing args" << std::endl;
        return 1;
    }
    
    // >> Mesh loading & preprocessing <<

    MyMesh m;

    int loadMask;
    if(ImporterOBJ<MyMesh>::Open(m, argv[1], loadMask) != ImporterOBJ<MyMesh>::E_NOERROR)
    {
        std::cout << "Error reading input file" << std::endl;
        return 1;
    }
    std::cout << "Loaded " << argv[1] << " with mask " << loadMask << std::endl;

    vcg::tri::RequireVFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
    vcg::tri::RequireFFAdjacency(m);
    vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
    vcg::tri::RequirePerVertexFlags<MyMesh>(m);
    vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceAdj(m);
    
    vcg::tri::RequirePerFaceNormal(m);

    // >> Gradient <<

    MyMesh::PerVertexAttributeHandle<vcg::Point3d> energyGrad = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<vcg::Point3d>(m, "Combinatorial energy gradient");
    double totEnergy;

    double dt;
    std::vector<double> dts;
    clock_t clockStart;

    // perform simple gradient descent with argv[2] steps and a stepsize of argv[3]
    int nSteps = atoi(argv[2]);
    float stepSize = atof(argv[3]);
    for(int step = 0; step < nSteps; step++)
    {
        clockStart = clock();
        vcg::tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);
        totEnergy = 0.0;

        // reset gradient of each vertex
        for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
            energyGrad[vIter].SetZero();

        combinatorialEnergyGrad(m, energyGrad, [&totEnergy](double energy, MyVertex* v) {
            totEnergy += energy;
        });

        for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
            (*vIter).P() -= (energyGrad[vIter] * stepSize);

        dt = ((double)(clock() - clockStart) / CLOCKS_PER_SEC);
        std::cout << "Step #" << step << ": Energy=" << totEnergy << "\tTime=" << dt << std::endl;
        dts.push_back(dt);    
    }

    dt = 0;
    for(double currDt : dts)
        dt += currDt;
    dt /= dts.size();
    std::cout << "Mean time per step: " << dt << std::endl;

    totEnergy = 0.0;
    for(MyMesh::VertexIterator vIter = m.vert.begin(); vIter != m.vert.end(); vIter++)
        totEnergy += localCombinatorialEnergy(&*vIter);

    std::cout << "Final energy: " << totEnergy << std::endl;

    vcg::tri::io::ExporterPLY<MyMesh>::Save(m, "out.ply", TriMask::IOM_NONE, false);

    return 0;
}