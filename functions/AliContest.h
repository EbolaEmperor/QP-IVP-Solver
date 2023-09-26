#include "IVP.h"
#include "functionFactory.h"

class AliContestFunc : public TimeFunction{
private:
    Real omega;
public:
    AliContestFunc(const Real &arg){
        omega = arg;
    }
    ColVector operator () (const ColVector &x, const Real &t) const{
        ColVector res(2);
        res(0) = x(1);
        res(1) = -(1 + 0.01 * cosq(omega*t)) * x(0);
        return res;
    }
};

static void registerAliContest(void)__attribute__((constructor));

void registerAliContest(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("Alibaba Contest", [](Real arg){ return (TimeFunction*) new AliContestFunc(arg); });
}