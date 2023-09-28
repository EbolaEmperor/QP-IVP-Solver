#include "IVP.h"
#include "functionFactory.h"

class nBodyFunc : public TimeFunction{
public:
    ColVector operator () (const ColVector &x, const Real &t) const{
        int n = x.size() / 4;
        ColVector res(x.size());
        for(int i = 0; i < n; i++){
            res(4*i) = x(4*i+2);
            res(4*i+1) = x(4*i+3);
            for(int j = 0; j < n; j++){
                if(i==j) continue;
                Real dx = x(4*i)-x(4*j), dy = x(4*i+1)-x(4*j+1);
                Real r3 = powq( sqr(dx) + sqr(dy), 1.5 );
                res(4*i+2) += -dx/r3;
                res(4*i+3) += -dy/r3;
            }
        }
        return res;
    }
};

static void registerNBody(void)__attribute__((constructor));

void registerNBody(){
    auto& factory = FunctionFactory::Instance();
    factory.registerFunction("N-Body Problem", [](Real arg){ return (TimeFunction*) new nBodyFunc(); });
}