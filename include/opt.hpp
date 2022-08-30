#ifndef OPT_HPP
#define OPT_HPP

#include "mesh.hpp"


/*
 * Gradient method optimizer interface
 */
class Optimizer
{
public:
    Optimizer(MyMesh& m, double stepSize);

    /*
     * Routine that must be called when the mesh geometry and/or topology changes
     */
    virtual void reset() = 0;

    /*
     * Perform an optimization step; returns true iff another step can be performed
     */
    virtual bool step() = 0;

    virtual void printStats() = 0;

    void updateGradientSqNorm();

    double getGradientSqNorm() { return gradSqNorm; }
    double getStepSize() { return stepSize; }
    double getEnergy() { return energy; }
    int getNFunEval() { return nFunEval; }

protected:
    MyMesh& m;
    AreaFaceAttrHandle fAttrArea;
    StarVertAttrHandle vAttrStar;
    GradientVertAttrHandle vAttrGrad;
    
    double stepSize;
    double gradSqNorm;
    double energy;
    int nFunEval;
};


/*
 * Gradient method optimization with fixed step size
 */
class FixedStepOpt : public Optimizer
{
public:
    FixedStepOpt(MyMesh& m,
                 int maxFunEval,
                 double eps,
                 double stepSize);

    void reset() override;
    bool step() override;
    void printStats() override;

private:
    int maxFunEval;
    double eps;
};


/*
 * Gradient method optimization with backtracking line search (Armijo condition)
 */
class BacktrackingOpt : public Optimizer
{
public:
    BacktrackingOpt(MyMesh& m,
                    int maxFunEval,
                    double eps,
                    double initialStepSize,
                    double tau,
                    double minStepSize = 1e-10,
                    double armijoM1 = 1e-4);

    void reset() override;
    bool step() override;
    void printStats() override;

private:
    std::vector<vcg::Point3d> tmpVP;
    int maxFunEval;
    double eps;
    double initialStepSize;
    double minStepSize;
    double tau;
    double armijoM1;
};


#endif