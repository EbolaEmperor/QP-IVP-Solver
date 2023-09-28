#include "IVP.h"
#include "functionFactory.h"

class MultiBodyInOneStarFunc : public TimeFunction{
public:
    ColVector operator () (const ColVector &x, const Real &t) const{
        int n = x.size() / 4;
        ColVector res(x.size());
        for(int i = 0; i < n; i++){
            res(4*i) = x(4*i+2);
            res(4*i+1) = x(4*i+3);
            Real r3 = powq( sqr(x(4*i)) + sqr(x(4*i+1)), 1.5 );
            res(4*i+2) = -x(4*i)/r3;
            res(4*i+3) = -x(4*i+1)/r3;
        }
        return res;
    }
};

static void registerMultiBodyInOneStar(void)__attribute__((constructor));

void registerMultiBodyInOneStar(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("N-Body 1-Star Problem", [](Real arg){ return (TimeFunction*) new MultiBodyInOneStarFunc(); });
}