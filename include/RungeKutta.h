#ifndef _RUNGE_KUTTA_
#define _RUNGE_KUTTA_

#include "IVP.h"

class ClassicalRKSolver : public TimeIntegrator_ConstStep{
protected:
    ColVector oneStepSolve(TimeFunction &f, const ColVector &x0, const ColVector &f0, const Real &t0, const Real &step);
    void solve(TimeFunction &f, const ColVector &x0, const Real &T, const int &gridSize);
};
static void registerClassicalRK(void)__attribute__((constructor));


class ImplicitRKSolver : public TimeIntegrator_ConstStep{
protected:
    virtual ColVector oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step) = 0;
    void solve(TimeFunction &f, const ColVector &x0, const Real &T, const int &gridSize);
};


class ESDIRKSolver : public ImplicitRKSolver{
public:
    ColVector oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step);
    ESDIRKSolver();
};
static void registerESDIRK(void)__attribute__((constructor));


class CollocationRKSolver_OneStepY{
protected:
    std::vector<ColVector> oneStepSolveY(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step, const Matrix &A, const ColVector &c);
};


class CollocationRKSolver : public ImplicitRKSolver{
protected:
    Matrix A;
    ColVector b, c;
    void computeAandB();
public:
    Matrix getA() const{return A;}
    ColVector getB() const{return b;}
    ColVector getC() const{return c;}
};


class ConstStepCollocationRKSolver : public CollocationRKSolver, public CollocationRKSolver_OneStepY{
protected:
    ColVector oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step);
};


class GaussLegendreRKSolver : public ConstStepCollocationRKSolver{
public:
    GaussLegendreRKSolver(const int &stage);
};
static void registerGaussLegendre(void)__attribute__((constructor));


class RadauIIARKSolver : public ConstStepCollocationRKSolver{
public:
    RadauIIARKSolver(const int &stage);
};
static void registerRadauIIA(void)__attribute__((constructor));


class SympleticRadauSolver : public ConstStepCollocationRKSolver{
public:
    SympleticRadauSolver();
};
static void registerSympleticRadau(void)__attribute__((constructor));


class AdaptiveRKSolver : public TimeIntegrator_VariativeStep{
protected:
    Matrix A;
    ColVector b1, c;
    virtual std::pair<ColVector, ColVector> oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step) = 0;
    void solve(TimeFunction &f, const ColVector &x0, const Real &T, const Real &eps);
};


class AdaptiveESDIRKSolver : public AdaptiveRKSolver{
private:
    ESDIRKSolver esdirk;
protected:
    std::pair<ColVector, ColVector> oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step);
public:
    AdaptiveESDIRKSolver();
};
static void registerAdaptiveESDIRK(void)__attribute__((constructor));

class EmbeddedESDIRKSolver : public AdaptiveRKSolver{
protected:
    Matrix A;
    ColVector b1, b2, c;
    std::pair<ColVector, ColVector> oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step);
public:
    EmbeddedESDIRKSolver(const int &p);
};
static void registerEmbeddedESDIRK(void)__attribute__((constructor));



class EmbeddedRKSolver : public AdaptiveRKSolver{
protected:
    ColVector b2;
    std::pair<ColVector, ColVector> oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step);
};


class FehlbergSolver : public EmbeddedRKSolver{
public:
    FehlbergSolver();
};
static void registerFehlberg(void)__attribute__((constructor));


class DormandPrinceSolver : public EmbeddedRKSolver{
public:
    DormandPrinceSolver();
};
static void registerDormandPrince(void)__attribute__((constructor));


class DormandPrince8Solver : public EmbeddedRKSolver{
public:
    DormandPrince8Solver();
};
static void registerDormandPrince8(void)__attribute__((constructor));


class AdaptiveCollocationRKSolver : public AdaptiveRKSolver, public CollocationRKSolver_OneStepY{
protected:
    ColVector trueOneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step);
    std::pair<ColVector, ColVector> oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step);
    void computeCoef(CollocationRKSolver *p);
};


class AdaptiveGaussLegendreRKSolver : public AdaptiveCollocationRKSolver{
public:
    AdaptiveGaussLegendreRKSolver(const int &stage);
};
static void registerAdaptiveGaussLegendre(void)__attribute__((constructor));


class AdaptiveRadauIIARKSolver : public AdaptiveCollocationRKSolver{
public:
    AdaptiveRadauIIARKSolver(const int &stage);
};
static void registerAdaptiveRadauIIA(void)__attribute__((constructor));


class AdaptiveSympleticRadauSolver : public AdaptiveCollocationRKSolver{
public:
    AdaptiveSympleticRadauSolver();
};
static void registerAdaptiveSympleticRadau(void)__attribute__((constructor));

#endif