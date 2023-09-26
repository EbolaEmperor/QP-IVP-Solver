#include "IVP.h"
#include "functionFactory.h"

class RestrictedThreeBodyFunc : public TimeFunction{
private:
    Real mu;
public:
    RestrictedThreeBodyFunc(const Real &_mu): mu(_mu){}
    ColVector operator () (const ColVector &x, const Real &t) const{
        ColVector res(6);
        res(0) = x(3);
        res(1) = x(4);
        res(2) = x(5);
        res(3) = 2*x(4) + x(0) - mu*(x(0)+mu-1) / powq(x(1)*x(1) + x(2)*x(2) + (x(0)+mu-1)*(x(0)+mu-1), 1.5) - (1-mu)*(x(0)+mu) / powq(x(1)*x(1) + x(2)*x(2) + (x(0)+mu)*(x(0)+mu), 1.5);
        res(4) = -2*x(3) + x(1) - mu*x(1) / powq(x(1)*x(1) + x(2)*x(2) + (x(0)+mu-1)*(x(0)+mu-1), 1.5) - (1-mu)*x(1) / powq(x(1)*x(1) + x(2)*x(2) + (x(0)+mu)*(x(0)+mu), 1.5);
        res(5) = -mu*x(2) / powq(x(1)*x(1) + x(2)*x(2) + (x(0)+mu-1)*(x(0)+mu-1), 1.5) - (1-mu)*x(2) / powq(x(1)*x(1) + x(2)*x(2) + (x(0)+mu)*(x(0)+mu), 1.5);
        return res;
    }
};

static void registerRestrictedThreeBody(void)__attribute__((constructor));

void registerRestrictedThreeBody(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Restricted 3-Body Problem", [](Real arg){ return (TimeFunction*) new RestrictedThreeBodyFunc(arg); });
}