#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
using namespace std;

long double h = .0001;
long double pi = 3.14159265359;
long double e = 2.71828182846;
long double E0 = .00001;


long double pFunction(long double yn,long double tn, long double p) {
    long double answer = 0;
    answer = pow(e,-1*p)-p-pow(yn,3)+(3*pow(e,-1*pow(tn,3)));
    return(answer);
}

long double pFunctionDeriv(long double yn, long double tn, long double p) {
    long double answer = 0;
    answer = (pFunction(yn,tn,(p+h))-pFunction(yn,tn,(p-h)))/(2*h);
    return(answer);
}

long double yPrime(long double yn, long double tn) {
    long double e = 1;
    long double pn = 0;
    while (e > E0) {
        pn = pn - (pFunction(yn,tn,pn))/(pFunctionDeriv(yn,tn,pn));
        e = pFunction(yn,tn,pn); 
    }

    return(pn);


}

void EulerEstimate(long double y0, long double x0, long double hE, long double n) {
    long double yn = y0;
    long double xn = 0;
    for(long double i = 1.0; i <= n; i = i+1) {
        xn = x0+(hE*i);
        yn = yn + (hE*yPrime(yn,xn));
        long double xIN = log2(xn);
        if(trunc(xIN)==xIN) {
            cout << "log2(tn) = " << xIN << "   y: " << yn << endl;
        }
        
    }
}

void HuenEstimate(long double y0, long double x0, long double hE, long double n) {
    long double yn = y0;
    long double xn = 0;
    for(long double i = 1.0; i <= n; i = i+1) {
        xn = x0+(hE*i);
        yn = yn + hE*(.25)*(yPrime(yn,xn) + (3*(yPrime( (xn+((2.0*hE)/3.0)),( yn+(2.0/3.0)*hE*yPrime(yn,xn)           )       ))));
        long double xIN = log2(xn);
        if(trunc(xIN)==xIN) {
            cout << "log2(tn) = " << xIN << "   y: " << yn << endl;
        }
        
    }
}








int main() {
    cout << "Eulers with h = 2^-10: " << endl;
    long double h1 = pow(2,-10);
    long double h2 = pow(2,-8);
    EulerEstimate(1,0,h1,16000);
    cout << "Huen with h = 2^-10: " << endl;
    HuenEstimate(1,0,h1,16000);
    cout << "Huen with h = 2^-8: " << endl;
    HuenEstimate(1,0,h2,16000);
    cout << "Huen with h = 2^-10: " << endl;
    EulerEstimate(1,0,h2,16000);
}