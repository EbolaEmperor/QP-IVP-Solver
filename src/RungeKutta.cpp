#include "RungeKutta.h"
#include "NonLinearSolver.h"
#include "Polynomial.h"

//-----------------------------------Classical RK Method-----------------------------------------

void registerClassicalRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Classical RK", [](int p){ return (TimeIntegrator*) new ClassicalRKSolver(); });
}

ColVector ClassicalRKSolver::oneStepSolve(TimeFunction &f, const ColVector &x0, const ColVector &f0, const Real &t0, const Real &step){
    ColVector RK_y1 = f0;
    ColVector RK_y2 = f(x0 + step/2*RK_y1, t0 + step/2);
    ColVector RK_y3 = f(x0 + step/2*RK_y2, t0 + step/2);
    ColVector RK_y4 = f(x0 + step*RK_y3, t0 + step);
    return x0 + step/6 * (RK_y1 + 2*RK_y2 + 2*RK_y3 + RK_y4);
}

void ClassicalRKSolver::solve(TimeFunction &f, const ColVector &x0, const Real &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    for(int i = 1; i <= n; i++){
        sol.push_back(oneStepSolve(f, sol[i-1], dsol[i-1], (i-1)*timeStep, timeStep));
        dsol.push_back(f(sol[i], i * timeStep));
    }
}

//-------------------------------------Implicit RK Method-----------------------------------------

void ImplicitRKSolver::solve(TimeFunction &f, const ColVector &x0, const Real &T, const int &n){
    maxTime = T;
    timeStep = maxTime/n;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    Real errC = -(Real)(1<<maxorder) / (1 - (1<<maxorder));
    for(int i = 1; i <= n; i++){
        sol.push_back(oneStepSolve(f, sol[i-1], (i-1)*timeStep, timeStep));
        dsol.push_back(f(sol[i],i*timeStep));
    }
}

//---------------------------------------ESDIRK Method-------------------------------------------

void registerESDIRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("ESDIRK", [](int p){ return (TimeIntegrator*) new ESDIRKSolver(); });
}

void registerAdaptiveESDIRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adaptive ESDIRK", [](int p){ return (TimeIntegrator*) new AdaptiveESDIRKSolver(); });
}

class ESDIRK_Function : public Function{
private:
    TimeFunction &f;
    ColVector x0;
    Real b, t;
public:
    ESDIRK_Function(TimeFunction &f, const ColVector &x0, const Real &b, const Real &t):
        f(f), x0(x0), b(b), t(t) {}
    ColVector operator () (const ColVector &x) const{
        return f(x0 + b * x, t) - x;
    }
};

ESDIRKSolver::ESDIRKSolver(){
    maxorder = 4;
}

ColVector ESDIRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step){
    static const Real aval[]={
        0,                              0,                          0,                          0,                      0,                  0,
        1.0Q/4.0Q,                        1.0Q/4.0Q,                    0,                          0,                      0,                  0,
        8611.0Q/62500.0Q,                 -1743.0Q/31250.0Q,            1.0Q/4.0Q,                    0,                      0,                  0,
        5012029.0Q/34652500.0Q,           -654441.0Q/2922500.0Q,        174375.0Q/388108.0Q,          1.0Q/4.0Q,                0,                  0,
        15267082809.0Q/155376265600.0Q,   -71443401.0Q/120774400.0Q,    730878875.0Q/902184768.0Q,    2285395.0Q/8070912.0Q,    1.0Q/4.0Q,            0,
        82889.0Q/524892.0Q,               0,                          15625.0Q/83664.0Q,            69875.0Q/102672.0Q,       -2260.0Q/8211.0Q,     1.0Q/4.0Q
    };
    static const Real cval[] = {
        0,      1.0Q/2.0Q,        83.0Q/250.0Q,     31.0Q/50.0Q,      17.0Q/20.0Q,      1.0Q
    };
    static const Matrix A(6, 6, aval);
    static const ColVector b(6, aval+30);
    static const ColVector c(6, cval);
    static const int s = 6;
    NonLinearSolver nlsolver;

    std::vector<ColVector> y;
    y.push_back(f(U0,t0));
    for(int i = 1; i < s; i++){
        ColVector constY = U0;
        for(int j = 0; j < i; j++)
            constY = constY + step * A(i,j) * y[j];
        Real t = t0 + c(i) * step;
        
        //针对线性函数的特殊优化（解热方程时有重要作用）
        if(f.isLinear()){
            ColVector _c(1);
            _c(0) = c(i);
            Matrix _a(1,1);
            _a(0,0) = A(i,i);
            y.push_back(f.solve(_a, _c, constY, t0, step));
            continue;
        }

        static const int MaxIter = 50;
        ColVector curY, nxtY = U0 + step * c(i) * y[0];
        int T = 0;
        do{
            curY = std::move(nxtY);
            nxtY = f(constY + step * A(i,i) * curY, t);
            if(++T==MaxIter) break;
            if(nxtY.maxnorm() > 1e100){
                T = MaxIter;
                break;
            }
        }while(relativeError(curY, nxtY) > 1e-30Q);
        
        if(T < MaxIter){
            y.push_back(nxtY);
        } else {
            // 通常情况下，不动点迭代是收敛的，如果不收敛，则需调用非线性方程组求解器
            ESDIRK_Function esf(f, constY, step*A(i,i), t);
            y.push_back(nlsolver.solve(esf, U0 + step*c(i)*y[0]));
        }
    }

    ColVector U1 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b(i) * y[i];
    }
    return U1;
}

AdaptiveESDIRKSolver::AdaptiveESDIRKSolver(){
    maxorder = 4;
}

std::pair<ColVector, ColVector> AdaptiveESDIRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step){
    ColVector U1 = esdirk.oneStepSolve(f, U0, t0, step);
    ColVector U2_mid = esdirk.oneStepSolve(f, U0, t0, step/2);
    ColVector U2 = esdirk.oneStepSolve(f, U2_mid, t0+step/2, step/2);
    return std::make_pair(U1,U2);
}

//--------------------------------------Embedded ESDIRK Methods----------------------------------------

void registerEmbeddedESDIRK(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Embedded ESDIRK", [](int p){ return (TimeIntegrator*) new EmbeddedESDIRKSolver(p); });
}

EmbeddedESDIRKSolver::EmbeddedESDIRKSolver(const int &p){
    if(p < 4 || p > 6){
        std::cerr << "The order of an Embedded ESDIRK method could only be 4 or 5." << std::endl;
        exit(-1);
    }
    if(p == 4){
        static const Real aval[]={
            0,                              0,                          0,                          0,                      0,                  0,                0,
            1.0Q/8.0Q,                     1.0Q/8.0Q,                 0,                          0,                      0,                  0,                 0,
            -39188347878.0Q/1513744654945.0Q,  -39188347878.0Q/1513744654945.0Q,  1.0Q/8.0Q,         0,                      0,                  0,                  0,
            1748874742213.0Q/5168247530883.0Q,  1748874742213.0Q/5168247530883.0Q,  -1748874742213.0Q/5795261096931.0Q,          1.0Q/8.0Q,                0,                  0,                 0,
            -6429340993097.0Q/17896796106705.0Q,  -6429340993097.0Q/17896796106705.0Q,  9711656375562.0Q/10370074603625.0Q,    1137589605079.0Q/3216875020685.0Q,    1.0Q/8.0Q,            0,                   0,
            405169606099.0Q/1734380148729.0Q,  405169606099.0Q/1734380148729.0Q,  -264468840649.0Q/6105657584947.0Q,     118647369377.0Q/6233854714037.0Q,       683008737625.0Q/4934655825458.0Q,     1.0Q/8.0Q,              0,
            -5649241495537.0Q/14093099002237.0Q,  -5649241495537.0Q/14093099002237.0Q,  5718691255176.0Q/6089204655961.0Q,  2199600963556.0Q/4241893152925.0Q,  8860614275765.0Q/11425531467341.0Q,  -3696041814078.0Q/6641566663007.0Q,       1.0Q/8.0Q
        };
        static const Real b2val[] = {
            -1517409284625.0Q/6267517876163.0Q,  -1517409284625.0Q/6267517876163.0Q,      8291371032348.0Q/12587291883523.0Q,        5328310281212.0Q/10646448185159.0Q,     5405006853541.0Q/7104492075037.0Q,      -4254786582061.0Q/7445269677723.0Q,      19.0Q/140.0Q
        };
        static const Real cval[] = {
            0,      1.0Q/4.0Q,        (2.0Q-sqrt(2.0Q))/8.0Q,     1.0Q/2.0Q,      395.0Q/567.0Q,      89.0Q/126.0Q,       1.0Q
        };
        A = Matrix(7, 7, aval);
        b1 = ColVector(7, aval+42);
        b2 = ColVector(7, b2val);
        c = ColVector(7, cval);
        maxorder = 4;
    } else if(p==5) {
        static const Real aval[]={
            0,                              0,                          0,                       0,                                                 0,                          0,                                    0,                              0,
            1.0Q/7.0Q,                     1.0Q/7.0Q,                       0,                       0,                                                 0,                          0,                                    0,                              0,
            1521428834970.0Q/8822750406821.0Q,  1521428834970.0Q/8822750406821.0Q,  1.0Q/7.0Q,         0,                                                 0,                          0,                                    0,                              0,
            5338711108027.0Q/29869763600956.0Q,  5338711108027.0Q/29869763600956.0Q,  1483184435021.0Q/6216373359362.0Q,          1.0Q/7.0Q,                0,                          0,                                    0,                              0,
            2264935805846.0Q/12599242299355.0Q,  2264935805846.0Q/12599242299355.0Q,  1330937762090.0Q/13140498839569.0Q,    -287786842865.0Q/17211061626069.0Q,    1.0Q/7.0Q,            0,                                    0,                              0,
            118352937080.0Q/527276862197.0Q,  118352937080.0Q/527276862197.0Q,  -2960446233093.0Q/7419588050389.0Q,     -3064256220847.0Q/46575910191280.0Q,       6010467311487.0Q/7886573591137.0Q,     1.0Q/7.0Q,              0,                              0,
            1134270183919.0Q/9703695183946.0Q,  1134270183919.0Q/9703695183946.0Q,  4862384331311.0Q/10104465681802.0Q,  1127469817207.0Q/2459314315538.0Q,  -9518066423555.0Q/11243131997224.0Q,  -811155580665.0Q/7490894181109.0Q,       1.0Q/7.0Q,              0,
            2162042939093.0Q/22873479087181.0Q,  2162042939093.0Q/22873479087181.0Q,  -4222515349147.0Q/9397994281350.0Q,  3431955516634.0Q/4748630552535.0Q,  -374165068070.0Q/9085231819471.0Q,  -1847934966618.0Q/8254951855109.0Q,  5186241678079.0Q/7861334770480.0Q,       1.0Q/7.0Q
        };
        static const Real b2val[] = {
            701879993119.0Q/7084679725724.0Q,  701879993119.0Q/7084679725724.0Q,  -8461269287478.0Q/14654112271769.0Q,      6612459227430.0Q/11388259134383.0Q,        2632441606103.0Q/12598871370240.0Q,     -2147694411931.0Q/10286892713802.0Q,      4103061625716.0Q/6371697724583.0Q,      36.0Q/233.0Q
        };
        static const Real cval[] = {
            0,      2.0Q/7.0Q,        (2.0Q+sqrt(2.0Q))/7.0Q,     150.0Q/203.0Q,      27.0Q/46.0Q,      473.0Q/532.0Q,       30.0Q/83.0Q,       1.0Q
        };
        A = Matrix(8, 8, aval);
        b1 = ColVector(8, aval+56);
        b2 = ColVector(8, b2val);
        c = ColVector(8, cval);
        maxorder = 5;
    } else {
        static const Real gma = 2.0Q/9.0Q;
        static const Real aval[]={
            0,                              0,                          0,                       0,                                                 0,                          0,                                    0,                              0, 0,
            gma,                            gma,                       0,                       0,                                                 0,                          0,                                    0,                              0, 0,
            1.0Q/9.0Q,                      -52295652026801.0Q/1014133226193379.0Q,  gma,         0,                                                 0,                          0,                                    0,                              0, 0,
            37633260247889.0Q/456511413219805.0Q,  -162541608159785.0Q/642690962402252.0Q,  186915148640310.0Q/408032288622937.0Q,          gma,                0,                          0,                                    0,                              0, 0,
            -37161579357179.0Q/532208945751958.0Q,  -211140841282847.0Q/266150973773621.0Q,  884359688045285.0Q/894827558443789.0Q,    845261567597837.0Q/1489150009616527.0Q,    gma,            0,                                    0,                              0, 0,
            32386175866773.0Q/281337331200713.0Q,  498042629717897.0Q/1553069719539220.0Q,  -73718535152787.0Q/262520491717733.0Q,     -147656452213061.0Q/931530156064788.0Q,       -16605385309793.0Q/2106054502776008.0Q,     gma,              0,                              0, 0,
            -38317091100349.0Q/1495803980405525.0Q,  233542892858682.0Q/880478953581929.0Q,  -281992829959331.0Q/709729395317651.0Q,  -52133614094227.0Q/895217507304839.0Q,  -9321507955616.0Q/673810579175161.0Q,  79481371174259.0Q/817241804646218.0Q,       gma,              0, 0,
            -486324380411713.0Q/1453057025607868.0Q,  -1085539098090580.0Q/1176943702490991.0Q,  370161554881539.0Q/461122320759884.0Q,  804017943088158.0Q/886363045286999.0Q,  -15204170533868.0Q/934878849212545.0Q,  -248215443403879.0Q/815097869999138.0Q,  339987959782520.0Q/552150039467091.0Q,       gma, 0,
            0, 0, 0, 281246836687281.0Q/672805784366875.0Q, 250674029546725.0Q/464056298040646.0Q, 88917245119922.0Q/798581755375683.0Q, 127306093275639.0Q/658941305589808.0Q, -319515475352107.0Q/658842144391777.0Q, gma
        };
        static const Real b2val[] = {
            -204006714482445.0Q/253120897457864.0Q,     0.0Q,     -818062434310719.0Q/743038324242217.0Q,  3176520686137389.0Q/1064235527052079.0Q,      -574817982095666.0Q/1374329821545869.0Q,        -507643245828272.0Q/1001056758847831.0Q,     2013538191006793.0Q/972919262949000.0Q,      352681731710820.0Q/726444701718347.0Q,      -12107714797721.0Q/746708658438760.0Q
        };
        static const Real cval[] = {
            0,      2.0Q*gma,        376327483029687.0Q/1335600577485745.0Q,     433625707911282.0Q/850513180247701.0Q,      183.0Q/200.0Q,      62409086037595.0Q/296036819031271.0Q,       81796628710131.0Q/911762868125288.0Q,     0.97Q,      1.0Q
        };
        A = Matrix(9, 9, aval);
        b1 = ColVector(9, aval+72);
        b2 = ColVector(9, b2val);
        c = ColVector(9, cval);
        maxorder = 6;
    }
}

std::pair<ColVector, ColVector> EmbeddedESDIRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step){
    NonLinearSolver nlsolver;
    std::vector<ColVector> y;
    const int s = c.size();
    y.push_back(f(U0,t0));

    for(int i = 1; i < s; i++){
        ColVector constY = U0;
        for(int j = 0; j < i; j++)
            constY = constY + step * A(i,j) * y[j];
        Real t = t0 + c(i) * step;
        
        //针对线性函数的特殊优化（解热方程时有重要作用）
        if(f.isLinear()){
            ColVector _c(1);
            _c(0) = c(i);
            Matrix _a(1,1);
            _a(0,0) = A(i,i);
            y.push_back(f.solve(_a, _c, constY, t0, step));
            continue;
        }

        static const int MaxIter = 100;
        ColVector curY, nxtY = U0 + step * c(i) * y[0];
        int T = 0;
        do{
            curY = std::move(nxtY);
            nxtY = f(constY + step * A(i,i) * curY, t);
            if(++T==MaxIter) break;
            if(nxtY.maxnorm() > 1e100){
                T = MaxIter;
                break;
            }
        }while(relativeError(curY, nxtY) > 1e-30Q);
        
        if(T < MaxIter){
            y.push_back(nxtY);
        } else {
            // 通常情况下，不动点迭代是收敛的，如果不收敛，则需调用非线性方程组求解器
            ESDIRK_Function esf(f, constY, step*A(i,i), t);
            y.push_back(nlsolver.solve(esf, U0 + step*c(i)*y[0]));
        }
    }

    ColVector U1 = U0, U2 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b1(i) * y[i];
        U2 = U2 + step * b2(i) * y[i];
    }
    return std::make_pair(U1,U2);
}

//--------------------------------------Gauss-Legendre RK Methods----------------------------------------

void registerGaussLegendre(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Gauss-Legendre", [](int p){ return (TimeIntegrator*) new GaussLegendreRKSolver(p); });
}

GaussLegendreRKSolver::GaussLegendreRKSolver(const int &s){
    if(s <= 0){
        std::cerr << "[Error] The stage of a Gauss-Legendre solver should be at least 1." << std::endl;
        exit(-1);
    }
    A = Matrix(s,s);
    b = ColVector(s);
    c = ColVector(s);
    
    Real *coef = new Real[s+1];
    for(int j = 0; j <= s; j++){
        coef[j] = sqr(fact(s))/fact(2*s) * binom(s,j) * binom(s+j,j);
        if((s-j)&1) coef[j] = -coef[j];
    }
    Polynomial poly(s, coef);
    auto roots = poly.roots();
    for(int i = 0; i < s; i++)
        c(i) = crealq(roots[i]);
    c.sort();
    computeAandB();
    maxorder = 2*s;
}

//--------------------------------------Radau IIA RK Methods----------------------------------------

void registerRadauIIA(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Radau-IIA", [](int p){ return (TimeIntegrator*) new RadauIIARKSolver(p); });
}

RadauIIARKSolver::RadauIIARKSolver(const int &s){
    if(s <= 0){
        std::cerr << "[Error] The stage of a Radau-IIA solver should be at least 1." << std::endl;
        exit(-1);
    }
    A = Matrix(s,s);
    b = ColVector(s);
    c = ColVector(s);

    Polynomial poly = pow(linearPolynomial(1,0),s-1) * pow(linearPolynomial(1,-1),s);
    for(int i = 0; i < s-1; i++) poly = poly.derivative();
    auto p = poly.roots();
    for(int i = 0; i < p.size(); i++)
        c(i) = crealq(p[i]);
    c.sort();
    computeAandB();
    maxorder = 2*s-1;
}

//--------------------------------------Sympletic Radau RK Method----------------------------------------

void registerSympleticRadau(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Sympletic Radau", [](int p){ return (TimeIntegrator*) new SympleticRadauSolver(); });
}

SympleticRadauSolver::SympleticRadauSolver(){
    const Real sq6 = sqrtq(6.0Q);
    const Real aval[] = {
        (16.0Q-sq6)/72.0Q,              (328.0Q-167.0Q*sq6)/1800.0Q,    (-2.0Q+3.0Q*sq6)/450.0Q,
        (328.0Q+167.0Q*sq6)/1800.0Q,    (16.0Q+sq6)/72.0Q,              (-2.0Q-3.0Q*sq6)/450.0Q,
        (85.0Q-10.0Q*sq6)/180.0Q,       (85.0Q+10.0Q*sq6)/180.0Q,       1.0Q/18.0Q
    };
    const Real bval[] = {
        (16.0Q-sq6)/36.0Q,  (16.0Q+sq6)/36.0Q,  1.0Q/9.0Q
    };
    const Real cval[] = {
        (4.0Q-sq6)/10.0Q,   (4.0Q+sq6)/10.0Q,   1.0Q
    };

    A = Matrix(3,3,aval);
    b = ColVector(3,bval);
    c = ColVector(3,cval);
    maxorder = 5;
}

//-------------------------Some Generic Functions for Adaptive RK Methods------------------------------

std::pair<ColVector, ColVector> EmbeddedRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step){
    static const int s = c.size();
    std::vector<ColVector> y;
    for(int i = 0; i < s; i++){
        ColVector x1 = U0;
        for(int j = 0; j < i; j++)
            x1 = x1 + step * A(i,j) * y[j];
        y.push_back(f(x1, t0 + c(i) * step));
    }

    ColVector U1 = U0, U2 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b1(i) * y[i];
        U2 = U2 + step * b2(i) * y[i];
    }
    return std::make_pair(U1, U2);
}

void AdaptiveRKSolver::solve(TimeFunction &f, const ColVector &x0, const Real &T, const Real &eps){
    static const Real rhomax = 5, rhomin = 0.2, rho = 0.2;
    maxTime = T;
    if(maxStep==-1) maxStep = T/100;
    int gridsize = 0;
    Real step = T/1e4Q, curtime = 0.0Q;
    sol.push_back(x0);
    dsol.push_back(f(x0,0));
    timePoint.push_back(0.0Q);
    while(curtime < T){
        if(curtime + step >= T) step = T - curtime;
        std::pair<ColVector, ColVector> pU = oneStepSolve(f, sol[gridsize], curtime, step);
        // std::cerr << (pU.first-pU.second).T() << " " << step << " " << curtime << std::endl;
        Real Eind = 0.0Q;
        for(int i = 0; i < pU.first.size(); i++){
            Eind += sqr( (pU.first(i) - pU.second(i)) / ( eps + fabsq(sol[gridsize](i)) * eps ) );
        }
        Eind = sqrtq(Eind / pU.first.size());
        if(Eind <= 1){
            curtime += step;
            sol.push_back(pU.first);
            dsol.push_back(f(pU.first, curtime));
            timePoint.push_back(curtime);
            gridsize++;
        }
        step = step * std::min( rhomax, std::max(rhomin, powq(rho/Eind, 1.0Q/(maxorder+1))) );
        if(step > maxStep) step = maxStep;
    }
}

//--------------------------------Fehlberg 4(5) Embedded RK Method-------------------------------------

void registerFehlberg(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Fehlberg", [](int p){ return (TimeIntegrator*) new FehlbergSolver(); });
}

FehlbergSolver::FehlbergSolver(){
    static const Real aval[] = {
        0,              0,              0,              0,              0,              0,
        1.0Q/4.0Q,        0,              0,              0,              0,              0,
        3.0Q/32.0Q,       9.0Q/32.0Q,       0,              0,              0,              0,
        1932.0Q/2197.0Q,  -7200.0Q/2197.0Q, 7296.0Q/2197.0Q, 0,              0,              0,
        439.0Q/216.0Q,    -8.0Q,           3680.0Q/513.0Q,   -845.0Q/4104.0Q,  0,              0,
        -8.0Q/27.0Q,      2.0Q,            -3544.0Q/2565.0Q, 1859.0Q/4104.0Q,  -11.0Q/40.0Q,     0,
    };
    static const Real b1val[] = {
        25.0Q/216.0Q, 0.0Q, 1408.0Q/2565.0Q, 2197.0Q/4104.0Q, -1.0Q/5.0Q, 0.0Q
    };
    static const Real b2val[] = {
        16.0Q/135.0Q, 0.0Q, 6656.0Q/12825.0Q, 28561.0Q/56430.0Q, -9.0Q/50.0Q, 2.0Q/55.0Q
    };
    static const Real cval[] = {
        0.0Q, 1.0Q/4.0Q, 3.0Q/8.0Q, 12.0Q/13.0Q, 1.0Q, 1.0Q/2.0Q
    };
    A = Matrix(6, 6, aval);
    b1 = ColVector(6, b1val);
    b2 = ColVector(6, b2val);
    c = ColVector(6, cval);
    maxorder = 4;
}

//--------------------------------Dormand-Prince 8(7) Embedded RK Method-------------------------------------

void registerDormandPrince8(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Dormand-Prince 8(7)", [](int p){ return (TimeIntegrator*) new DormandPrince8Solver(); });
}

DormandPrince8Solver::DormandPrince8Solver(){
    static const Real aval[] = {
        0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        1.0Q/18.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        1.0Q/48.0Q, 1.0Q/16.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        1.0Q/32.0Q, 0.0Q, 3.0Q/32.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        5.0Q/16.0Q, 0.0Q, -75.0Q/64.0Q, 75.0Q/64.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        3.0Q/80.0Q, 0.0Q, 0.0Q, 3.0Q/16.0Q, 3.0Q/20.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        29443841.0Q/614563906.0Q, 0.0Q, 0.0Q, 77736538.0Q/692538347.0Q, -28693883.0Q/1125000000.0Q, 23124283.0Q/1800000000.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        16016141.0Q/946692911.0Q, 0.0Q, 0.0Q, 61564180.0Q/158732637.0Q, 22789713.0Q/633445777.0Q, 545815736.0Q/2771057229.0Q, -180193667.0Q/1043307555.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        39632708.0Q/573591083.0Q, 0.0Q, 0.0Q, -433636366.0Q/683701615.0Q, -421739975.0Q/2616292301.0Q, 100302831.0Q/723423059.0Q, 790204164.0Q/839813087.0Q, 800635310.0Q/3783071287.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        246121993.0Q/1340847787.0Q, 0.0Q, 0.0Q, -37695042795.0Q/15268766246.0Q, -309121744.0Q/1061227803.0Q, -12992083.0Q/490766935.0Q, 6005943493.0Q/2108947869.0Q, 393006217.0Q/1396673457.0Q, 123872331.0Q/1001029789.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q,
        -1028468189.0Q/846180014.0Q, 0.0Q, 0.0Q, 8478235783.0Q/508512852.0Q, 1311729495.0Q/1432422823.0Q, -10304129995.0Q/1701304382.0Q, -48777925059.0Q/3047939560.0Q, 15336726248.0Q/1032824649.0Q, -45442868181.0Q/3398467696.0Q, 3065993473.0Q/597172653.0Q, 0.0Q, 0.0Q, 0.0Q,
        185892177.0Q/718116043.0Q, 0.0Q, 0.0Q, -3185094517.0Q/667107341.0Q, -477755414.0Q/1098053517.0Q, -703635378.0Q/230739211.0Q, 5731566787.0Q/1027545527.0Q, 5232866602.0Q/850066563.0Q, -4093664535.0Q/808688257.0Q, 3962137247.0Q/1805957418.0Q, 65686358.0Q/487910083.0Q, 0.0Q,0.0Q,
        403863854.0Q/491063109.0Q, 0.0Q, 0.0Q, -5068492393.0Q/434740067.0Q, -411421997.0Q/543043805.0Q, 652783627.0Q/914296604.0Q, 11173962825.0Q/925320556.0Q, -13158990841.0Q/6184727034.0Q, 3936647629.0Q/1978049680.0Q, -160528059.0Q/685178525.0Q, 248638103.0Q/1413531060.0Q, 0.0Q, 0.0Q,
    };
    static const Real b1val[] = {
        14005451.0Q/335480064.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, -59238493.0Q/1068277825.0Q, 181606767.0Q/758867731.0Q, 561292985.0Q/797845732.0Q, -1041891430.0Q/1371343529.0Q, 760417239.0Q/1151165299.0Q, 118820643.0Q/751138087.0Q, -528747749.0Q/2220607170.0Q, 1.0Q/4.0Q
    };
    static const Real b2val[] = {
        13451932.0Q/455176623.0Q, 0.0Q, 0.0Q, 0.0Q, 0.0Q, -808719846.0Q/976000145.0Q, 1757004468.0Q/5645159321.0Q, 656045339.0Q/265891186.0Q, -3867574721.0Q/1518517206.0Q, 465885868.0Q/322736535.0Q, 53011238.0Q/667516719.0Q, 2.0Q/45.0Q, 0.0Q
    };
    static const Real cval[] = {
        0.0Q, 1.0Q/18.0Q, 1.0Q/12.0Q, 1.0Q/8.0Q, 5.0Q/16.0Q, 3.0Q/8.0Q, 59.0Q/400.0Q, 93.0Q/200.0Q, 5490023248.0Q/9719169821.0Q, 13.0Q/20.0Q, 1201146811.0Q/1299019798.0Q, 1.0Q, 1.0Q
    };
    A = Matrix(13, 13, aval);
    b1 = ColVector(13, b1val);
    b2 = ColVector(13, b2val);
    c = ColVector(13, cval);
    maxorder = 7;
}

//--------------------------------Dormand-Prince 5(4) Embedded RK Method-------------------------------------

void registerDormandPrince(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Dormand-Prince", [](int p){ return (TimeIntegrator*) new DormandPrinceSolver(); });
}

DormandPrinceSolver::DormandPrinceSolver(){
    static const Real aval[] = {
        0,              0,                  0,                0,              0,                   0,           0,
        1.0Q/5.0Q,        0,                  0,                0,              0,                   0,           0,
        3.0Q/40.0Q,       9.0Q/40.0Q,           0,                0,              0,                   0,           0,
        44.0Q/45.0Q,      -56.0Q/15.0Q,         32.0Q/9.0Q,         0,              0,                   0,           0,
        19372.0Q/6561.0Q, -25360.0Q/2187.0Q,    64448.0Q/6561.0Q,   -212.0Q/729.0Q,   0,                   0,           0,
        9017.0Q/3168.0Q,  -355.0Q/33.0Q,        46732.0Q/5247.0Q,   49.0Q/176.0Q,     -5103.0Q/18656.0Q,     0,           0,
        35.0Q/384.0Q,     0.0Q,                500.0Q/1113.0Q,     125.0Q/192.0Q,    -2187.0Q/6784.0Q,      11.0Q/84.0Q,   0
    };
    static const Real b1val[] = {
        35.0Q/384.0Q,     0.0Q,                500.0Q/1113.0Q,     125.0Q/192.0Q,    -2187.0Q/6784.0Q,      11.0Q/84.0Q,   0
    };
    static const Real b2val[] = {
        5179.0Q/57600.0Q,     0.0Q,    7571.0Q/16695.0Q,     393.0Q/640.0Q,    -92097.0Q/339200.0Q,  187.0Q/2100.0Q,   1.0Q/40.0Q
    };
    static const Real cval[] = {
        0.0Q, 1.0Q/5.0Q, 3.0Q/10.0Q, 4.0Q/5.0Q, 8.0Q/9.0Q, 1.0Q, 1.0Q
    };
    A = Matrix(7, 7, aval);
    b1 = ColVector(7, b1val);
    b2 = ColVector(7, b2val);
    c = ColVector(7, cval);
    maxorder = 4;
}

//-------------------------------------Adaptive Gauss-Legendre RK Method--------------------------------------------

void registerAdaptiveGaussLegendre(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adaptive Gauss-Legendre", [](int p){ return (TimeIntegrator*) new AdaptiveGaussLegendreRKSolver(p); });
}

AdaptiveGaussLegendreRKSolver::AdaptiveGaussLegendreRKSolver(const int &s){
    maxorder = 2*s;
    CollocationRKSolver *p = new GaussLegendreRKSolver(s);
    computeCoef(p);
}

//-------------------------------------Adaptive Radau IIA RK Method--------------------------------------------

void registerAdaptiveRadauIIA(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adaptive Radau-IIA", [](int p){ return (TimeIntegrator*) new AdaptiveRadauIIARKSolver(p); });
}

AdaptiveRadauIIARKSolver::AdaptiveRadauIIARKSolver(const int &s){
    maxorder = 2*s - 1;
    CollocationRKSolver *p = new RadauIIARKSolver(s);
    computeCoef(p);
}

//-------------------------------------Adaptive Sympletic Radau RK Method--------------------------------------------

void registerAdaptiveSympleticRadau(){
    auto& factory = TimeIntegratorFactory::Instance();
    factory.registerTimeIntegrator("Adaptive Sympletic Radau", [](int p){ return (TimeIntegrator*) new AdaptiveSympleticRadauSolver(); });
}

AdaptiveSympleticRadauSolver::AdaptiveSympleticRadauSolver(){
    maxorder = 5;
    CollocationRKSolver *p = new SympleticRadauSolver();
    computeCoef(p);
}

//-------------------------Some Generic Functions for Collocation Methods---------------------------

void CollocationRKSolver::computeAandB(){
    const int s = c.size();
    std::vector<Polynomial> l;
    for(int i = 0; i < s; i++){
        Polynomial p = constPolynomial(1);
        for(int j = 0; j < s; j++){
            if(i == j) continue;
            p = p * linearPolynomial(1, -c(j)) / (c(i) - c(j));
        }
        l.push_back(p.integral());
    }

    for(int i = 0; i < s; i++){
        b(i) = l[i](1) - l[i](0);
        for(int j = 0; j < s; j++)
            A(i,j) = l[j](c(i)) - l[j](0);
    }
}

class CollocationSolver_Function : public Function{
private:
    TimeFunction &f;
    const Matrix &A;
    const ColVector &c;
    const ColVector &U0;
    const Real step, t0;
public:
    CollocationSolver_Function(TimeFunction &_f, const Matrix &_A,  const ColVector &_c, const ColVector &_U0, const Real &_step, const Real &_t0):
        f(_f), A(_A), c(_c), U0(_U0), step(_step), t0(_t0){}
    ColVector operator () (const ColVector &x) const{
        const int s = c.size(), m = x.size()/s;
        ColVector res(s*m);
        std::vector<ColVector> y;
        for(int i = 0; i < s; i++)
            y.push_back(x.getSubmatrix(i*m,i*m+m-1,0,0));
        for(int i = 0; i < s; i++){
            ColVector Xk = U0;
            for(int j = 0; j < s; j++)
                Xk = Xk + step * A(i,j) * y[j];
            res.setSubmatrix(i*m, 0, f(Xk, t0 + c(i)*step) - y[i]);
        }
        return res;
    }
};

std::vector<ColVector> CollocationRKSolver_OneStepY::oneStepSolveY(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step, const Matrix &A, const ColVector &c){
    NonLinearSolver nlsolver;
    std::vector<ColVector> y;
    const int s = c.size(), m = U0.size();
    y.push_back(f(U0,t0));
    for(int i = 1; i < s; i++)
        y.push_back(y[0]);
    int T = 0;
    static const int MaxIter = 100;
    while(!f.isLinear() && ++T < MaxIter){
        Real maxerr = 0;
        for(int i = 0; i < s; i++){
            ColVector sumY = U0;
            for(int j = 0; j < s; j++)
                sumY = sumY + step * A(i,j) * y[j];
            ColVector nxtY = f(sumY, t0 + c(i) * step);
            maxerr = std::max(maxerr, relativeError(nxtY,y[i]));
            y[i] = nxtY;
            if(y[i].maxnorm() > 1e100){
                T = MaxIter-1;
                break;
            }
        }
        if(maxerr < 1e-30Q) break;
    }
    if(f.isLinear() || T == MaxIter){
        // 通常情况下，不动点迭代是收敛的，如果不收敛，则需调用非线性方程组求解器
        ColVector resY;
        if(!f.isLinear()){
            CollocationSolver_Function glf(f, A, c, U0, step, t0);
            ColVector initY(s*m), y0 = f(U0,t0);
            for(int i = 0; i < s; i++)
                initY.setSubmatrix(i*m, 0, y0);
            resY = nlsolver.solve(glf, initY);
        } else {
            // 对于线性函数f，可以直接求解线性方程组，不需要迭代，这将大大优化求解速度
            resY = f.solve(A, c, U0, t0, step);
        }
        for(int i = 0; i < s; i++)
            y[i] = resY.getSubmatrix(i*m, i*m+m-1, 0, 0);
    }
    return y;
}

ColVector ConstStepCollocationRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step){
    const int s = c.size();
    std::vector<ColVector> y = oneStepSolveY(f, U0, t0, step, A, c);
    ColVector U1 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b(i) * y[i];
    }
    return U1;
}

ColVector AdaptiveCollocationRKSolver::trueOneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step){
    const int s = c.size();
    std::vector<ColVector> y = oneStepSolveY(f, U0, t0, step, A, c);
    ColVector U1 = U0;
    for(int i = 0; i < s; i++){
        U1 = U1 + step * b1(i) * y[i];
    }
    return U1;
}

std::pair<ColVector, ColVector> AdaptiveCollocationRKSolver::oneStepSolve(TimeFunction &f, const ColVector &U0, const Real &t0, const Real &step){
    ColVector U1 = trueOneStepSolve(f, U0, t0, step);
    ColVector U2_mid = trueOneStepSolve(f, U0, t0, step/2);
    ColVector U2 = trueOneStepSolve(f, U2_mid, t0+step/2, step/2);
    return std::make_pair(U1,U2);
}

void AdaptiveCollocationRKSolver::computeCoef(CollocationRKSolver *p){
    A = p->getA();
    b1 = p->getB();
    c = p->getC();
}