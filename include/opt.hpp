#ifndef OPT_HPP
#define OPT_HPP

#include "mesh_matrix.hpp"


class Optimizer
{
public:
    Optimizer(int nVertices, double stepSize) :
        G(nVertices, 3),
        stepSize(stepSize),
        gradNorm(-1),
        currEnergy(-1),
        nFunEval(0)
    {}

    virtual bool step(Eigen::Ref<Matrix3Xd> V,
                      const Eigen::Ref<const Matrix3Xi>& F,
                      Eigen::Ref<Matrix3Xd> N,
                      Eigen::Ref<ArrayXd> A,
                      const Eigen::Ref<const MatrixXi>& S,
                      const Eigen::Ref<const ArrayXb>& B) = 0;

    virtual void printStats() = 0;

    double getGradientNorm() { return gradNorm; }
    double getCurrStepSize() { return stepSize; }
    double getCurrEnergy() { return currEnergy; }
    int getNFunEval() { return nFunEval; }

protected:
    Matrix3Xd G;
    double stepSize;
    double gradNorm;
    double currEnergy;
    int nFunEval;
};


// Gradient method optimization with fixed step size
class FixedStepOpt : public Optimizer
{
public:
    FixedStepOpt(int nVertices,
                 int maxFunEval,
                 double eps,
                 double stepSize);

    bool step(Eigen::Ref<Matrix3Xd> V,
              const Eigen::Ref<const Matrix3Xi>& F,
              Eigen::Ref<Matrix3Xd> N,
              Eigen::Ref<ArrayXd> A,
              const Eigen::Ref<const MatrixXi>& S,
              const Eigen::Ref<const ArrayXb>& B) override;

    void printStats() override;

private:
    int maxFunEval;
    double eps;
};

// Gradient method optimization with backtracking line search (Armijo condition)
class BacktrackingOpt : public Optimizer
{
public:
    BacktrackingOpt(int nVertices,
                    int maxFunEval,
                    double eps,
                    double initialStepSize = 1.0,
                    double minStepSize = 1e-16,
                    double tau = 0.9,
                    double armijoM1 = 0.0001);

    bool step(Eigen::Ref<Matrix3Xd> V,
            const Eigen::Ref<const Matrix3Xi>& F,
            Eigen::Ref<Matrix3Xd> N,
            Eigen::Ref<ArrayXd> A,
            const Eigen::Ref<const MatrixXi>& S,
            const Eigen::Ref<const ArrayXb>& B) override;

    void printStats() override;

private:
    int maxFunEval;
    double eps;
    double initialStepSize;
    double minStepSize;
    double tau;
    double armijoM1;

    Matrix3Xd tmpV;
};

// Gradient method optimization with Lewis and Overton line search (Armijo + strong Wolfe conditions)
class LewisOvertonOpt : public Optimizer
{
public:
    LewisOvertonOpt(int nVertices,
                    int maxFunEval,
                    double eps,
                    double initialStepSize = 1.0,
                    double minStepSize = 1e-12,
                    double armijoM1 = 0.0001,
                    double wolfeM3 = 0.99);

    bool step(Eigen::Ref<Matrix3Xd> V,
            const Eigen::Ref<const Matrix3Xi>& F,
            Eigen::Ref<Matrix3Xd> N,
            Eigen::Ref<ArrayXd> A,
            const Eigen::Ref<const MatrixXi>& S,
            const Eigen::Ref<const ArrayXb>& B) override;

    void printStats() override;

private:
    int maxFunEval;
    double eps;
    double initialStepSize;
    double minStepSize;
    double armijoM1;
    double wolfeM3;

    Matrix3Xd tmpV;
    Matrix3Xd tmpG;
};


#endif