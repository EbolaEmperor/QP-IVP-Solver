#include "IVP.h"
#include "functionFactory.h"

class VanDerPolFunc : public TimeFunction{
public:
    ColVector operator () (const ColVector &x, const Real &t) const{
        ColVector res(2);
        res(0) = x(1);
        res(1) = 1000*(1-x(0)*x(0))*x(1) - x(0);
        return res;
    }
};

static void registerVanDerPol(void)__attribute__((constructor));

void registerVanDerPol(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Van der Pol Problem", [](Real arg){ return (TimeFunction*) new VanDerPolFunc(); });
}