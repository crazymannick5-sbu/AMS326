#include <iostream>
#include <cmath>

int main() {

}

double fFunction(double x){
    double A = 0;
    double B = 0;
    B=-1*pow(x,2);
    A = exp(B);

}

double BisectionMain(double a1, double b1) {
    double x1 = (a1 + b1)/2.0;
    fa1 = fFunction(a1);
    fx1 = fFunction(x1);
    Mult = fa1 * fx1;
    if(abs(a1-b1) < .00001) {return((a1+b1)/2)    } else {
        if(Mult < 0) {
            return(BisectionMain(a1, x1));
        } else {
            return(BisectionMain(x1,b1));
        }
    }

}

