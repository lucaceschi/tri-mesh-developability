#ifndef OPT_HPP
#define OPT_HPP

#include "mesh.hpp"


class Optimizer
{
public:
    Optimizer(MyMesh& m,
              double stepSize) :
        m(m),
        stepSize(stepSize),
        nFunEval(0)
    {}

    virtual void reset() = 0;
    virtual bool step() = 0;
    virtual void printStats() = 0;

    void updateGradientSqNorm(GradientVertAttrHandle vAttrGrad);
    double getGradientSqNorm() { return gradSqNorm; }
    double getStepSize() { return stepSize; }
    double getEnergy() { return energy; }
    int getNFunEval() { return nFunEval; }

protected:
    MyMesh& m;
    double stepSize;
    double gradSqNorm;
    double energy;
    int nFunEval;
};


// Gradient method optimization with fixed step size
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
    AreaFaceAttrHandle fAttrArea;
    StarVertAttrHandle vAttrStar;
    GradientVertAttrHandle vAttrGrad;
    int maxFunEval;
    double eps;
};


// Gradient method optimization with backtracking line search (Armijo condition)
class BacktrackingOpt : public Optimizer
{
public:
    BacktrackingOpt(MyMesh& m,
                    int maxFunEval,
                    double eps,
                    double initialStepSize = 0.8,
                    double minStepSize = 1e-10,
                    double tau = 0.8,
                    double armijoM1 = 1e-4);

    void reset() override;
    bool step() override;
    void printStats() override;

private:
    MyMesh tmpMesh;
    AreaFaceAttrHandle fAttrArea;
    StarVertAttrHandle vAttrStar;
    GradientVertAttrHandle vAttrGrad;
    int maxFunEval;
    double eps;
    double initialStepSize;
    double minStepSize;
    double tau;
    double armijoM1;
};


// Gradient method optimization with Lewis and Overton line search (Armijo + strong Wolfe conditions)
/*class LewisOvertonOpt : public Optimizer
{
public:
    LewisOvertonOpt(MyMesh& m,
                    int maxFunEval,
                    double eps,
                    double initialStepSize = 0.8,
                    double minStepSize = 1e-10,
                    double armijoM1 = 0.0001,
                    double wolfeM3 = 0.99);

    void reset() override;
    bool step() override;
    void printStats() override;

private:
    int maxFunEval;
    double eps;
    double initialStepSize;
    double minStepSize;
    double armijoM1;
    double wolfeM3;
    MyMesh tmpMesh;
    GradientVertAttrHandle tmpvAttrGrad;
};*/


#endif