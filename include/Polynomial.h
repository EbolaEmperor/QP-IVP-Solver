#ifndef _POLYNOMIAL_H_
#define _POLYNOMIAL_H_

#include "matrix.h"

// 一个不完整的多项式运算库

class Polynomial{
private:
    int n;
    Real *a;
    void removeLeadingZero();
public:
    Polynomial(): n(-1), a(nullptr){}
    Polynomial(const int &_n);
    Polynomial(const int &_n, const Real *p);
    Polynomial(const Polynomial &rhs);
    Polynomial& operator = (const Polynomial &rhs);
    Polynomial operator + (const Polynomial &rhs) const;
    Polynomial operator * (const Polynomial &rhs) const;
    friend Polynomial pow(Polynomial lhs, int b);
    Polynomial operator / (const Real &c) const;
    Polynomial operator / (const Polynomial &rhs) const;
    Polynomial derivative() const;
    Polynomial integral() const;
    std::vector<Complex> roots();
    Real operator () (const Real &x) const;
    Real& coef(const int &i);
    Real coef(const int &i) const;
    void print() const;
    friend std::ostream& operator << (std::ostream &out, const Polynomial &p);
};

Polynomial constPolynomial(const Real &a);
Polynomial linearPolynomial(const Real &a, const Real &b);
Polynomial HermiteInterpolation32(const Real &x0, const Real &x1, const Real &y0, const Real &dy0, const Real &y1, const Real &dy1);

#endif