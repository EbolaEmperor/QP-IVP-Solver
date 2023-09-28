#include "IVP.h"
#include "functionFactory.h"

class DoublePendulumFunc : public TimeFunction{
private:
    Real m2;
public:
    DoublePendulumFunc(const Real &x): m2(x) {}
    ColVector operator () (const ColVector &x, const Real &t) const{
        ColVector res(4);
        res(0) = x(2);
        res(1) = x(3);
        Matrix A(2,2);
        A(0,0) = 1+m2;
        A(0,1) = m2*0.8*cos(x(0)-x(1));
        A(1,0) = 0.8*cos(x(0)-x(1));
        A(1,1) = 0.8;
        ColVector rhs(2);
        rhs(0) = (1+m2)*9.8*cos(x(0)) - m2*0.8*sin(x(0)-x(1))*x(3)*x(3);
        rhs(1) = sin(x(0)-x(1))*x(2)*x(2) + 9.8*cos(x(1));
        ColVector sol = A.solve(rhs);
        res(2) = sol(0);
        res(3) = sol(1);
        return res;
    }
};

static void registerDoublePendulum(void)__attribute__((constructor));

void registerDoublePendulum(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Double Pendulum Problem", [](Real arg){ return (TimeFunction*) new DoublePendulumFunc(arg); });
}