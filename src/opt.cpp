#include "opt.hpp"

#include "energy_grad.hpp"
#include <iostream>
#include <iomanip>


FixedStepOpt::FixedStepOpt(int nVertices,
                           int maxFunEval,
                           double eps,
                           double stepSize) :
    Optimizer(nVertices, stepSize),
    maxFunEval(maxFunEval),
    eps(eps)
{}

bool FixedStepOpt::step(Eigen::Ref<Matrix3Xd> V,
                        const Eigen::Ref<const Matrix3Xi>& F,
                        Eigen::Ref<Matrix3Xd> N,
                        Eigen::Ref<ArrayXd> A,
                        const Eigen::Ref<const MatrixXi>& S,
                        const Eigen::Ref<const ArrayXb>& B)
{
    if(nFunEval >= maxFunEval || (nFunEval > 0 && gradNorm <= eps))
        return false;
    
    computeNormals(V, F, N, A);

    double totEnergy = 0.0;

    combinatorialEnergyGrad(V, F, N, A, S, B, G, [&totEnergy](double localEnergy, size_t v) {
        totEnergy += localEnergy;
    });
    nFunEval++;

    currEnergy = totEnergy;
    gradNorm = G.norm();

    V -= (G * stepSize);

    std::cout << "[FixedStepOpt] nFunEval=" << getNFunEval()
        << std::setw(20) << "stepSize=" << getCurrStepSize()
        << std::setw(20) << "gradNorm=" << getGradientNorm()
        << std::setw(20) << "energy=" << getCurrEnergy()
        << std::endl;

    return true;
}

void FixedStepOpt::printStats()
{
    std::cout << "[FixedStepOpt] nFunEval=" << getNFunEval()
        << std::setw(20) << "stepSize=" << getCurrStepSize()
        << std::setw(20) << "gradNorm=" << getGradientNorm()
        << std::setw(20) << "energy=" << getCurrEnergy()
        << std::endl;
}


BacktrackingOpt::BacktrackingOpt(int nVertices,
                                 int maxFunEval,
                                 double eps,
                                 double initialStepSize,
                                 double minStepSize,
                                 double tau,
                                 double armijoM1) :
    Optimizer(nVertices, initialStepSize),
    maxFunEval(maxFunEval),
    eps(eps),
    initialStepSize(initialStepSize),
    minStepSize(minStepSize),
    tau(tau),
    armijoM1(armijoM1),
    tmpV(nVertices, 3)
{}

bool BacktrackingOpt::step(Eigen::Ref<Matrix3Xd> V,
                           const Eigen::Ref<const Matrix3Xi>& F,
                           Eigen::Ref<Matrix3Xd> N,
                           Eigen::Ref<ArrayXd> A,
                           const Eigen::Ref<const MatrixXi>& S,
                           const Eigen::Ref<const ArrayXb>& B)
{    
    if(nFunEval >= maxFunEval || (nFunEval > 0 && gradNorm <= eps))
        return false;

    if(gradNorm == -1)
    {
        // first step
        double totEnergy = 0.0;
        computeNormals(V, F, N, A);
        combinatorialEnergyGrad(V, F, N, A, S, B, G, [&totEnergy](double localEnergy, size_t v) {
            totEnergy += localEnergy;
        });
        currEnergy = totEnergy;
    }

    gradNorm = G.norm();
    stepSize = initialStepSize;
    
    double lineSearchEnergy;
    double currTomographyDeriv = -(G.array().pow(2)).sum(); // == dot prod between G and search direction -G

    while(stepSize > minStepSize)
    {
        tmpV = V - G * stepSize;

        computeNormals(tmpV, F, N, A);
        lineSearchEnergy = combinatorialEnergy(N, S, B);
        nFunEval++;

        // check Armijo condition
        if(lineSearchEnergy <= currEnergy + armijoM1 * stepSize * currTomographyDeriv)
            break;

        if(nFunEval >= maxFunEval)
            return false;

        stepSize *= tau;
    }

    if(stepSize < minStepSize)
        stepSize = minStepSize;

    V = tmpV;
    currEnergy = lineSearchEnergy;
    combinatorialEnergyGrad(V, F, N, A, S, B, G);
    nFunEval++;

    return true;
}

void BacktrackingOpt::printStats()
{
    std::cout << "[BacktrackingOpt] nFunEval=" << getNFunEval()
        << std::setw(20) << "stepSize=" << getCurrStepSize()
        << std::setw(20) << "gradNorm=" << getGradientNorm()
        << std::setw(20) << "energy=" << getCurrEnergy()
        << std::endl;
}


LewisOvertonOpt::LewisOvertonOpt(int nVertices,
                                 int maxFunEval,
                                 double eps,
                                 double initialStepSize,
                                 double minStepSize,
                                 double armijoM1,
                                 double wolfeM3) :
    Optimizer(nVertices, initialStepSize),
    maxFunEval(maxFunEval),
    eps(eps),
    initialStepSize(initialStepSize),
    minStepSize(minStepSize),
    armijoM1(armijoM1),
    wolfeM3(wolfeM3),
    tmpV(nVertices, 3),
    tmpG(nVertices, 3)
{}

bool LewisOvertonOpt::step(Eigen::Ref<Matrix3Xd> V,
                           const Eigen::Ref<const Matrix3Xi>& F,
                           Eigen::Ref<Matrix3Xd> N,
                           Eigen::Ref<ArrayXd> A,
                           const Eigen::Ref<const MatrixXi>& S,
                           const Eigen::Ref<const ArrayXb>& B)
{
    if(nFunEval >= maxFunEval || (nFunEval > 0 && gradNorm <= eps))
        return false;

    if(gradNorm == -1)
    {
        // first step
        double totEnergy = 0.0;
        computeNormals(V, F, N, A);
        combinatorialEnergyGrad(V, F, N, A, S, B, G, [&totEnergy](double localEnergy, size_t v) {
            totEnergy += localEnergy;
        });
        currEnergy = totEnergy;
    }

    double lineSearchEnergy;
    double currTomographyDeriv = -(G.array().pow(2)).sum(); // == dot prod between G and search dir -G
    double tmpTomographyDeriv;
    double stepSizeLowerBound = 0;
    double stepSizeUpperBound = std::numeric_limits<double>::max();

    stepSize = initialStepSize;
    while(stepSize > minStepSize)
    {
        tmpV = V - G * stepSize;

        computeNormals(tmpV, F, N, A);
        lineSearchEnergy = 0;
        combinatorialEnergyGrad(tmpV, F, N, A, S, B, tmpG, [&lineSearchEnergy](double localEnergy, size_t v) {
            lineSearchEnergy += localEnergy;
        });
        nFunEval++;

        // check if Armijo condition fails
        if(lineSearchEnergy > currEnergy + armijoM1 * stepSize * currTomographyDeriv)
            stepSizeUpperBound = stepSize;
        else
        {
            tmpTomographyDeriv = -(G.array() * tmpG.array()).sum(); // == dot prod between tmpG and search dir -G
            // check if strong Wolfe condition fails
            if(tmpTomographyDeriv < wolfeM3 * currTomographyDeriv)
                stepSizeLowerBound = stepSize;
            else
                break;
        }

        if(stepSizeUpperBound != std::numeric_limits<double>::max())
            stepSize = (stepSizeLowerBound + stepSizeUpperBound) / 2.0;
        else
            stepSize = 2 * stepSizeLowerBound;

        if(nFunEval >= maxFunEval)
            return false;
    }

    if(stepSize < minStepSize)
        stepSize = minStepSize;

    V -= (G * stepSize);
    lineSearchEnergy = 0;
    combinatorialEnergyGrad(V, F, N, A, S, B, G, [&lineSearchEnergy](double localEnergy, size_t v) {
        lineSearchEnergy += localEnergy;
    });
    nFunEval++;
    currEnergy = lineSearchEnergy;
    gradNorm = G.norm();

    return true;
}

void LewisOvertonOpt::printStats()
{
    std::cout << "[LewisOvertonOpt] nFunEval=" << getNFunEval()
        << std::setw(20) << "stepSize=" << getCurrStepSize()
        << std::setw(20) << "gradNorm=" << getGradientNorm()
        << std::setw(20) << "energy=" << getCurrEnergy()
        << std::endl;
}