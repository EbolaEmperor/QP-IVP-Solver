#include "Polynomial.h"

Polynomial::Polynomial(const int &_n){
    n = _n;
    a = new Real[n+1];
    memset(a, 0, sizeof(Real)*(n+1));
}

Polynomial::Polynomial(const int &_n, const Real *p){
    n = _n;
    a = new Real[n+1];
    memcpy(a, p, sizeof(Real)*(n+1));
}

Polynomial::Polynomial(const Polynomial & rhs){
    n = rhs.n;
    a = new Real[n+1];
    memcpy(a, rhs.a, sizeof(Real)*(n+1));
}

Polynomial& Polynomial::operator= (const Polynomial &rhs){
    n = rhs.n;
    a = new Real[n+1];
    memcpy(a, rhs.a, sizeof(Real)*(n+1));
    return *this;
}

Polynomial constPolynomial(const Real &a){
    Real p[] = {a};
    return Polynomial(0, p);
}

Polynomial linearPolynomial(const Real &a, const Real &b){
    Real p[] = {b, a};
    return Polynomial(1, p);
}

Polynomial Polynomial::operator + (const Polynomial &rhs) const{
    Polynomial res(std::max(n,rhs.n));
    for(int i = 0; i <= n; i++)
        res.coef(i) += coef(i);
    for(int i = 0; i <= rhs.n; i++)
        res.coef(i) += rhs.coef(i);
    return res;
}

Polynomial Polynomial::operator * (const Polynomial &rhs) const{
    Polynomial res(n+rhs.n);
    for(int i = 0; i <= n; i++)
        for(int j = 0; j <= rhs.n; j++)
            res.a[i+j] += a[i] * rhs.a[j];
    return res;
}

Polynomial pow(Polynomial a, int b){
    Polynomial res(1);
    res.coef(0) = 1;
    for(;b; b >>= 1, a = a * a)
        if(b & 1) res = res * a;
    return res;
}

Polynomial Polynomial::operator / (const Real &c) const{
    Polynomial res = (*this);
    for(int i = 0; i <= n; i++)
        res.a[i] /= c;
    return res;
}

Polynomial Polynomial::operator / (const Polynomial &rhs) const{
    Polynomial res(n-rhs.n), lhs(*this);
    for(int i = n; i >= rhs.n; i--){
        Real tmp = lhs.coef(i) / rhs.coef(rhs.n);
        res.coef(i-rhs.n) = tmp;
        for(int j = 0; j <= rhs.n; j++)
            lhs.coef(i-j) -= rhs.coef(rhs.n-j) * tmp;
    }
    return res;
}

Real Polynomial::operator() (const Real &x) const{
    Real rx = 1.0, res = 0.0;
    for(int i = 0; i <= n; i++){
        res += a[i] * rx;
        rx *= x;
    }
    return res;
}

Real& Polynomial::coef(const int &i){
    return a[i];
}

Real Polynomial::coef(const int &i) const{
    return a[i];
}

Polynomial Polynomial::derivative() const{
    Polynomial res(n-1);
    for(int i = 0; i < n; i++)
        res.a[i] = a[i+1] * (i+1);
    return res;
}

Polynomial Polynomial::integral() const{
    Polynomial res(n+1);
    for(int i = 0; i <= n; i++)
        res.a[i+1] = a[i]/(i+1);
    return res;
}

std::vector<Complex> Polynomial::roots(){
    removeLeadingZero();
    Matrix A(n,n);
    for(int i = 1; i < n; i++)
        A(i,i-1) = 1;
    for(int i = 0; i < n; i++)
        A(i,n-1) = -a[i]/a[n];
    return A.eigen();
}

void Polynomial::removeLeadingZero(){
    while(n>=0 && a[n]==0) n--;
    Real *b = a;
    a = n>=0 ? new Real[n+1] : nullptr;
    if(n>=0) memcpy(a, b, sizeof(Real)*(n+1));
    delete [] b;
}

void Polynomial::print() const{
    for(int i = n; i >= 0; i--)
        std::cout << a[i] << " ";
    std::cout << std::endl;
}

Polynomial HermiteInterpolation32(const Real &x0, const Real &x1, const Real &y0, const Real &dy0, const Real &y1, const Real &dy1){
    Matrix A(4,4);
    ColVector b(4);
    A(0,0) = 1; A(0,1) = x0; A(0,2) = x0*x0; A(0,3) = x0*x0*x0;
    A(1,0) = 0; A(1,1) = 1;  A(1,2) = 2*x0;  A(1,3) = 3*x0*x0;
    A(2,0) = 1; A(2,1) = x1; A(2,2) = x1*x1; A(2,3) = x1*x1*x1;
    A(3,0) = 0; A(3,1) = 1;  A(3,2) = 2*x1;  A(3,3) = 3*x1*x1;
    b(0) = y0;
    b(1) = dy0;
    b(2) = y1;
    b(3) = dy1;
    ColVector coef = A.solve(b);
    Real c[] = {coef(0), coef(1), coef(2), coef(3)};
    return Polynomial(3, c);
}

std::ostream& operator << (std::ostream &out, const Polynomial &p){
    int sti = p.n;
    while(sti >= 0 && p.a[sti] == 0) sti--;
    for(int i = sti; i >= 0; i--){
        if(p.a[i]==0) continue;
        if(i<sti && p.a[i]>=0) out << '+';
        if(p.a[i]!=1 || i==0){
            out << p.a[i];
            if(i) out << '*';
        }
        if(i>=2) out << "x^" << i;
        if(i==1) out << "x";
    }
    if(sti < 0) out << '0';
    return out;
}
