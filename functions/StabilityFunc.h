#include "IVP.h"
#include "functionFactory.h"

class StabilityFunc : public TimeFunction{
private:
    static Real lambda;
public:
    ColVector trueSol(const Real &t) const{
        ColVector res(1);
        res(0) = expq(lambda*t);
        return res;
    }
    ColVector operator () (const ColVector &x, const Real &t) const{
        ColVector res(1);
        res(0) = lambda * x(0);
        return res;
    }
};

Real StabilityFunc::lambda = -1e6;

static void registerStability(void)__attribute__((constructor));

void registerStability(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Stability Problem", [](Real arg){ return (TimeFunction*) new StabilityFunc(); });
}