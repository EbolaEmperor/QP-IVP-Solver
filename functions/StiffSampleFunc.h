#include "IVP.h"
#include "functionFactory.h"

class StiffSampleFunc : public TimeFunction{
private:
    static Real lambda, eta;
public:
    ColVector trueSol (const ColVector &x, const Real &t) const{
        ColVector res(1);
        res(0) = expq(lambda * t) * (eta - 1) + cosq(t);
        return res;
    }
    ColVector operator () (const ColVector &x, const Real &t) const{
        ColVector res(1);
        res(0) = lambda * ( x(0) - cosq(t) ) - sinq(t);
        return res;
    }
};

Real StiffSampleFunc::lambda = -1e6;
Real StiffSampleFunc::eta = 1.5;

static void registerStiffSample(void)__attribute__((constructor));

void registerStiffSample(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Stiff Sample Problem", [](Real arg){ return (TimeFunction*) new StiffSampleFunc(); });
}