#include "opt.hpp"

#include "energy_grad.hpp"
#include <iostream>
#include <iomanip>


/*/ OPTIMIZER INTERFACE /*/

Optimizer::Optimizer(MyMesh& m, double stepSize) :
    m(m),
    stepSize(stepSize),
    nFunEval(0)
{
    vAttrStar = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<Star>(m, std::string("Star"));
    vAttrGrad = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<vcg::Point3d>(m, std::string("Gradient"));
    fAttrArea = vcg::tri::Allocator<MyMesh>::GetPerFaceAttribute<double>(m, std::string("Area"));
}


void Optimizer::updateGradientSqNorm()
{
    gradSqNorm = 0.0;
    
    for(size_t v = 0; v < m.VN(); v++)
        for(int i = 0; i < 3; i++)
            gradSqNorm += std::pow(vAttrGrad[v][i], 2);
}


/*/ FIXED STEP SIZE OPTIMIZER /*/

FixedStepOpt::FixedStepOpt(MyMesh& m,
                           int maxFunEval,
                           double eps,
                           double stepSize) :
    Optimizer(m, stepSize),
    maxFunEval(maxFunEval),
    eps(eps)
{
    reset();
}


void FixedStepOpt::reset()
{
    updateFaceStars(m, vAttrStar);
    updateNormalsAndAreas(m, fAttrArea);
    energy = combinatorialEnergyGrad(m, fAttrArea, vAttrStar, vAttrGrad);
    updateGradientSqNorm();
}


bool FixedStepOpt::step()
{
    if(nFunEval >= maxFunEval || gradSqNorm <= eps)
        return false;
    
    for(size_t v = 0; v < m.VN(); v++)
        m.vert[v].P() -= (vAttrGrad[v] * stepSize);

    updateNormalsAndAreas(m, fAttrArea);
    energy = combinatorialEnergyGrad(m, fAttrArea, vAttrStar, vAttrGrad);
    updateGradientSqNorm();
    nFunEval++;

    return true;
}


void FixedStepOpt::printStats()
{
    std::cout << "[FixedStepOpt] nFunEval=" << getNFunEval()
        << std::setw(20) << "stepSize=" << getStepSize()
        << std::setw(20) << "gradSqNorm=" << getGradientSqNorm()
        << std::setw(20) << "energy=" << getEnergy()
        << std::endl;
}


/*/ BACKTRACKING OPTIMIZER /*/

BacktrackingOpt::BacktrackingOpt(MyMesh& m,
                                 int maxFunEval,
                                 double eps,
                                 double initialStepSize,
                                 double tau,
                                 double minStepSize,
                                 double armijoM1) :
    Optimizer(m, initialStepSize),
    maxFunEval(maxFunEval),
    eps(eps),
    initialStepSize(initialStepSize),
    tau(tau),
    minStepSize(minStepSize),
    armijoM1(armijoM1)
{
    reset();
}


void BacktrackingOpt::reset()
{
    tmpVP.clear();
    tmpVP.reserve(m.vert.size());
    for(size_t v = 0; v < m.vert.size(); v++)
        tmpVP.push_back(m.vert[v].cP());
        
    updateFaceStars(m, vAttrStar);
    updateNormalsAndAreas(m, fAttrArea);
    energy = combinatorialEnergyGrad(m, fAttrArea, vAttrStar, vAttrGrad);
    updateGradientSqNorm();
}


bool BacktrackingOpt::step()
{
    if(nFunEval >= maxFunEval || gradSqNorm <= eps)
        return false;
    
    double LS_energy, LS_stepSize;
    // compute the current tomography deriv as the dot prod between grad and search direction -grad
    // == minus the squared l2 norm of the gradient
    double tomographyDeriv = -gradSqNorm;

    for(LS_stepSize = initialStepSize; LS_stepSize > minStepSize; LS_stepSize *= tau)
    {
        for(size_t v = 0; v < m.vert.size(); v++)
            m.vert[v].P() = tmpVP[v] - vAttrGrad[v] * LS_stepSize;

        updateNormalsAndAreas(m, fAttrArea);
        LS_energy = combinatorialEnergy(m, vAttrStar);
        nFunEval++;

        // check Armijo condition
        if(LS_energy <= energy + armijoM1 * LS_stepSize * tomographyDeriv)
            break;

        if(nFunEval >= maxFunEval)
        {
            for(size_t v = 0; v < m.vert.size(); v++)
                m.vert[v].P() = tmpVP[v];
            return false;
        }
    }

    for(size_t v = 0; v < m.vert.size(); v++)
        tmpVP[v] = m.vert[v].cP();

    stepSize = LS_stepSize;
    energy = LS_energy;
    combinatorialEnergyGrad(m, fAttrArea, vAttrStar, vAttrGrad);
    updateGradientSqNorm();
    nFunEval++;

    return true;
}


void BacktrackingOpt::printStats()
{
    std::cout << "[BacktrackingOpt] nFunEval=" << getNFunEval()
        << std::setw(20) << "stepSize=" << getStepSize()
        << std::setw(20) << "gradNorm=" << getGradientSqNorm()
        << std::setw(20) << "energy=" << getEnergy()
        << std::endl;
}