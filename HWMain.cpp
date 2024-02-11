#include <iostream>
#include <cmath>

int main() {
    double answer = 0;
    answer = BisectionMain(0,1);
}

/*double fFunction(double x){
    double A = 0;
    double B = 0;
    B=-1*pow(x,2.0);
    A = exp(B);
    return(A-pow(x,3.0))

}*/

double fFunction(double x) {
    return(pow(x,3.0)+x-1)
}

double BisectionMain(double a1, double b1) {
    double x1 = (a1 + b1)/2.0;
    double fa1 = fFunction(a1);
    double fx1 = fFunction(x1);
    double Mult = fa1 * fx1;
    double ans = 0;
    
    if(abs(a1-b1) < .00001) {ans = (a1+b1)/2; std::cout << ans << std::endl;    } else {
        if(Mult < 0) {
            ans = BisectionMain(a1, x1);
        } else {
            ans = BisectionMain(x1,b1);
        }
    }
    return(ans);

}

