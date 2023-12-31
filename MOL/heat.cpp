#include <bits/stdc++.h>
#include "IVP.h"
#include "json.h"
using namespace std;

const Real PI = acosq(-1);

class HeatEquation : public TimeFunction{
private:
    Real h;
    // 在这里设置边值条件
    Real b0(const Real &t) const{
        return 0.0;
    }
    Real b1(const Real &t) const{
        return 0.0;
    }
public:
    HeatEquation(const int &n){
        __isLinear = true;
        h = 1.0 / n;
    }
    ColVector operator () (const ColVector &U, const Real &t) const{
        ColVector U1(U.size());
        for(int i = 0; i < U.size(); i++){
            U1(i) = -2*U(i)/(h*h);
            if(i-1 >= 0) U1(i) += U(i-1)/(h*h);
            if(i+1 < U.size()) U1(i) += U(i+1)/(h*h);
            if(i==0) U1(i) += b0(t) / (h*h);
            if(i+1==U.size()) U1(i) += b1(t) / (h*h);
        }
        return U1;
    }
    ColVector solve(const Matrix &a, const ColVector &c, const ColVector &U0, const Real &t0, const Real &k) const {
        int s = c.size(), m = U0.size();
        ColVector rhs(s*m);
        for(int i = 0; i < s; i++)
            rhs.setSubmatrix(i*m, 0, (*this)(U0, t0+c(i)*k));
        Matrix coef(s*m, s*m);
        for(int i = 0; i < s; i++)
            for(int j = 0; j < s; j++)
                for(int l = 0; l < m; l++){
                    coef(i*m+l, j*m+l) = a(i,j)*k*2.0/(h*h);
                    if(i==j) coef(i*m+l, j*m+l) += 1;
                    if(l) coef(i*m+l, j*m+l-1) = -a(i,j)*k*1.0/(h*h);
                    if(l<m-1) coef(i*m+l, j*m+l+1) = -a(i,j)*k*1.0/(h*h);
                }
        return coef.solve(rhs);
    }
};

// 在这里设置初值条件
Real g0(const Real &x){
    if(0.45 <= x && x < 0.5) return 20*(x-0.45);
    else if(0.5 <= x && x < 0.55) return -20*(x-0.55);
    else return 0;
}

// Fourier Solver不支持非零边值条件
class FourierSolver{
private:
    int maxn;
    vector<Real> A;

    Real F(const Real &x) const{
        return 2 * initial(x) * sinq(maxn*PI*x);
    }

    Real squareIntegrate(const Real &l, const Real &r) const{
        Real mid = (l + r) / 2;
        return (F(l) + 4*F(mid) + F(r)) * (r - l) / 6;
    }

    // Simpson自适应积分
    Real integrate(const Real &l, const Real &r, const Real &V, const Real &eps) const{
        Real mid = (l + r) / 2;
        Real L = squareIntegrate(l,mid);
        Real R = squareIntegrate(mid,r);
        if(fabsq(L+R-V) <= 15*eps) return L + R + (L+R-V) / 15.0;
        return integrate(l,mid,L,eps) + integrate(mid,r,R,eps);
    }

    Real integrate(const Real &l, const Real &r, const Real &eps) const{
        // 将积分区间随机划分为9个小区间，防止出现因为5个等分点为0导致整个积分返回0的情况
        int points[10], m;
        for(int i = 1; i < 9; i++) points[i] = rand();
        points[0] = 0; points[9] = RAND_MAX;
        sort(points, points+10);
        m = unique(points, points+10) - points;
        Real sum = 0;
        for(int i = 0; i < m-1; i++){
            Real hl = l + (Real)points[i]/RAND_MAX * (r-l);
            Real hr = l + (Real)points[i+1]/RAND_MAX * (r-l);
            sum += integrate(hl, hr, squareIntegrate(hl,hr), max(1e-4Q*eps, 1e-30Q));
        }
        return sum;
    }

public:
    FourierSolver(const int &N){
        cout << "--------------------------------------------------------------------------------" << endl;
        cout << "Solving with Fourier solver..." << endl;
        maxn = 0;
        do{
            ++maxn;
            A.push_back(integrate(0, 1, 1e-15));
        } while(maxn < N);
        cout << maxn << " terms has been retained in the Fourier series." << endl;
    }

    Real initial(const Real &x) const{
        return g0(x);
    }

    Real operator () (const Real &x, const Real &t) const{
        Real sum = 0;
        for(int k = 1; k <= maxn; k++){
            sum += A[k-1] * expq(-k*k*PI*PI*t) * sinq(k*PI*x);
        }
        return sum;
    }

    void output(const string &fname, const Real &T, const Real &h){
        cout << "--------------------------------------------------------------------------------" << endl;
        ofstream fout(fname);
        fout << setprecision(12) << T << " ";
        for(Real x = h; x <= 1-h+1e-29Q; x += h){
            fout << (*this)(x, T) << " ";
        }
        fout << endl;
        fout.close();
        cout << "Output: Results has been saved to " << fname << endl;
    }

    void denseDiscreteOutput(const string &fname, const Real &step, const Real &T, const Real &h){
        cout << "--------------------------------------------------------------------------------" << endl;
        ofstream fout(fname);
        fout << setprecision(12);
        for(Real x = 0; x <= 1+1e-29Q; x += h){
            fout << initial(x) << " ";
        }
        fout << endl;
        for(Real t = step; t <= T+1e-29Q; t += step){
            fout << t << " ";
            for(Real x = h; x <= 1-h+1e-29Q; x += h){
                fout << (*this)(x, T) << " ";
            }
            fout << endl;
        }
        fout.close();
        cout << "Dense-Discrete Output: Results has been saved to " << fname << endl;
    }
};

int main(int argc, char* argv[]){
    Json::Reader reader;
    Json::Value problem;
    ifstream ifs;
    if(argc < 2){
        cerr << "Please provide input file name!" << endl;
        exit(-1);
    }
    ifs.open(argv[1]);
    if(!ifs.is_open()){
        cerr << "cannot read file " << argv[1] << endl;
        exit(-1);
    }
    if(!reader.parse(ifs, problem)){
        cerr << "parse error" << endl;
        exit(-1);
    }

    int n = problem["Space Section"].asInt();
    Real h = 1.0 / n;
    Real T = problem["End Time"].asDouble();

    if(problem["Method"].asString() == "Fourier"){
        FourierSolver solver(problem["Fourier Terms"].asInt());
        if(problem["Output"].asBool())
            solver.output("result.txt", T, h);
        if(problem["Dense-Discrete Output"].asBool())
            solver.denseDiscreteOutput("result-dense.txt", problem["Dense-Discrete Output Step"].asDouble(), T, h);
        return 0;
    }

    ColVector U0(n-1);
    for(int i = 1; i < n; i++)
        U0(i-1) = g0(i*h);
    HeatEquation f(n);

    auto& factory = TimeIntegratorFactory::Instance();
    auto solver = factory.createTimeIntegrator(problem["Method"].asString(), problem["Order"].asInt());
    if(problem.isMember("Time Section")){
        solver->solveWithInfo(f, U0, T, problem["Time Section"].asInt());
    } else {
        if(problem.isMember("Max Step"))
            solver->setMaxStep(problem["Max Step"].asDouble());
        solver->solveWithInfo(f, U0, T, (Real)problem["Tolerance"].asDouble());
    }
    
    if(problem["Output"].asBool())
        solver->output("result.txt");
    if(problem["Dense-Discrete Output"].asBool())
        solver->denseDiscreteOutput("result-dense.txt", problem["Dense-Discrete Output Step"].asDouble());
    return 0;
}