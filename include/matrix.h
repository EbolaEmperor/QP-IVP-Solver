/***************************************************************
 *
 * 这是一个矩阵运算库，为了方便以后设计算法更加简洁，特编写以用
 * 版本号：v1.1.0
 * 
 * copyright © 2023 Wenchong Huang, All rights reserved.
 *
 **************************************************************/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
#include <complex>
#include <quadmath.h>

typedef __float128 Real;
typedef __complex128 Complex;

std::ostream& operator << (std::ostream& out, const Real &a);

const Real mat_eps = 1000 * FLT128_EPSILON;

class Matrix;
class ColVector;
class RowVector;

class Matrix{
protected:
    Real *a;
public:
    int n, m;
    Matrix();
    Matrix(const int &_n);
    Matrix(const int &_n, const int &_m);
    Matrix(const Matrix &A);
    Matrix(const Real *p, const int &_n);
    Matrix(const int &_n, const int &_m, const Real *p);
    ~Matrix();
    bool empty() const;

    Matrix & operator = (const Matrix & rhs);
    Matrix & operator = (Matrix && rhs);

    const Real operator () (const int &r, const int &c) const;
    Real & operator () (const int &r, const int &c);
    const Real element(const int &r, const int &c) const;
    Real & element(const int &r, const int &c);

    friend Matrix diag(const Matrix &A);

    void setSubmatrix(const int &r, const int &c, const Matrix &rhs);
    Matrix getSubmatrix(const int &r1, const int &r2, const int &c1, const int &c2) const;
    Matrix reshape(const int &_n, const int &_m) const;

    RowVector getRow(const int &r) const;
    ColVector getCol(const int &c) const;

    Matrix operator + (const Real &x) const;
    Matrix operator - (const Real &x) const;
    Matrix operator + (const Matrix &B) const;
    Matrix operator - () const;
    Matrix operator - (const Matrix &B) const;
    Matrix operator * (const Matrix &B) const;
    Matrix operator / (const Real &p) const;
    Matrix T() const;

    Real vecnorm(const Real &p) const;
    Real maxnorm() const;
    void swaprow(const int &r1, const int &r2);
    void swapcol(const int &r1, const int &r2);

    friend Matrix solve(Matrix A, Matrix b);
    ColVector solve(const ColVector &b) const;
    Real det() const;
    Matrix inv() const;
    Matrix rref() const;
    void FGdecompose(Matrix &F, Matrix &G) const;
    Matrix pinv() const;
    Real sqrsum() const;

//下面是1.1.0版本新增函数，用于求解矩阵的特征值，预计下一版本添加反幂法求特征向量
public:
    std::vector<Complex> eigen() const;
private:
    Matrix realSchur() const;
    std::pair<Matrix,Matrix> hessenberg() const;
    std::pair<ColVector,Real> householder() const;
    std::pair<Matrix,Matrix> RealQR() const;
    std::pair<Matrix,Matrix> getQR() const;
    bool isComplexEigen() const;
    std::pair<Complex,Complex> getComplexEigen() const;
};

class RowVector: public Matrix{
public:
    RowVector(): Matrix() {};
    RowVector(const int &n): Matrix(1,n) {};
    RowVector(const int &n, const Real *p): Matrix(1,n,p) {};
    RowVector(const Matrix &rhs);
    int size() const;
    const Real operator ()(const int &x) const;
    Real & operator () (const int &x);
    RowVector operator + (const RowVector &rhs) const;
    RowVector operator - (const RowVector &rhs) const;
    RowVector operator - () const;
    ColVector T() const;
};

class ColVector: public Matrix{
public:
    ColVector(): Matrix() {};
    ColVector(const int &n): Matrix(n,1) {};
    ColVector(const int &n, const Real *p): Matrix(n,1,p) {};
    ColVector(const Matrix &rhs);
    int size() const;
    void sort();
    const Real operator ()(const int &x) const;
    Real & operator () (const int &x);
    ColVector operator + (const ColVector &rhs) const;
    ColVector operator - (const ColVector &rhs) const;
    ColVector operator - () const;
    RowVector T() const;
};

Matrix hilbert(const int &n);
Matrix zeros(const int &n, const int &m);
Matrix ones(const int &n, const int &m);
Matrix eye(const int &n);
Real value(const Matrix &A);
ColVector zeroCol(const int &n);
RowVector zeroRow(const int &n);
int sgn(const Real &x);

//----------------------Matrix相关函数---------------------------
Matrix operator * (const Real &k, const Matrix &x);
Matrix abs(const Matrix &A);
Matrix log2(const Matrix &A);
Real max(const Matrix &A);
Real sum(const Matrix &A);
std::ostream& operator << (std::ostream& out, const Matrix &A);
Matrix dotdiv(const Matrix &a, const Matrix &b);
Matrix solveLowerTriangular(const Matrix &A, const Matrix &b);
Matrix solveUpperTriangular(const Matrix &A, const Matrix &b);
ColVector CG_solve(const Matrix &A, const ColVector &b);
ColVector CG_solve(const Matrix &A, const ColVector &b, const Real err);
ColVector CG_solve(const Matrix &A, const ColVector &b, const Real err, ColVector x);
Real det(const Matrix &A);
Matrix inv(const Matrix &A);
Matrix choleskyImproved(const Matrix &A);
Matrix solveByLDL(const Matrix &A, const Matrix &b);
Matrix gillMurray(Matrix A);
Matrix solveByLDL_GM(const Matrix &A, const Matrix &b);
Matrix pinv(const Matrix &A);
Matrix mergeCol(const Matrix &A, const Matrix &B);
Matrix mergeRow(const Matrix &A, const Matrix &B);
Matrix min(const Matrix &A, const Matrix &B);
Matrix max(const Matrix &A, const Matrix &B);
Real vecnorm(const Matrix &A, const Real &p);
Real vecnorm(const Matrix &A);
Real maxnorm(const Matrix &A);
Real relativeError(const Matrix &A, const Matrix &B);
Matrix randMatrix(const int &n, const int &m);
Matrix randInvertibleMatrix(const int &n);

//----------------------Row/ColVector相关函数----------------------
RowVector operator * (const Real &k, const RowVector &x);
ColVector operator * (const Real &k, const ColVector &x);
ColVector operator * (const Matrix &A, const ColVector &x);
RowVector operator * (const RowVector &x, const Matrix &A);
Real operator * (const RowVector &r, const ColVector &c);

//----------------------基本算术扩展----------------------------
Real fact(const int &n);
Real sqr(const Real &x);
Real binom(const int &n, const int &k);

#endif