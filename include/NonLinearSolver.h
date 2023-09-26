#ifndef _NON_LIN_SOLVER_H_
#define _NON_LIN_SOLVER_H_

#include "matrix.h"

class Function{
public:
    virtual ColVector operator () (const ColVector &x) const = 0;
    Matrix jacobi(const ColVector &x) const;
    Real fval(const ColVector &x) const;
    ColVector gradval(const ColVector &x) const;
};


class NonLinearSolver{
private:
    // Armijo非精确搜索
    Real search(Function &f, const ColVector &x0, const ColVector &direct, const ColVector &grad);
public:
    // 非线性方程求解器（牛顿法）
    ColVector solve(Function &F, ColVector current);
    ColVector solve(Function &F, ColVector current, const Real &err);
    ColVector solve(Function &F, ColVector current, const Real &err, const int &MAXN);
};

#endif