#include "IVP.h"
#include "functionFactory.h"

class HamiltonianFunc : public TimeFunction{
public:
    ColVector operator () (const ColVector &x, const Real &t) const{
        ColVector res(4);
        const Real r = x(2)*x(2) + x(3)*x(3);
        res(0) = - x(2) / powq(r, 1.5Q) - 0.015Q * x(2) / powq(r, 2.5Q);
        res(1) = - x(3) / powq(r, 1.5Q) - 0.015Q * x(3) / powq(r, 2.5Q);
        res(2) = x(0);
        res(3) = x(1);
        return res;
    }
};

static void registerHamiltonian(void)__attribute__((constructor));

void registerHamiltonian(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Hamiltonian", [](Real arg){ return (TimeFunction*) new HamiltonianFunc(); });
}