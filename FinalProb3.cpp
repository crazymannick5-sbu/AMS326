#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
#include <complex>
#include <iomanip>
using namespace std; 
double e = 2.71828182;
double pi = 3.141592659;
double y0og = 12;
double h= .0001;
double kNaught = 3.5;

void displayLoadingBar(double value) {
    const int barWidth = 70;

    // Calculate progress ratio
    double ratio = std::abs(value) / kNaught;
    int progress = static_cast<int>(ratio * barWidth);

    // Output the loading bar
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < progress) std::cout << "=";
        else std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(2) << (ratio * 100.0) << "%\r";
    std::cout.flush();
}
double Yprime(double y, double k) {
    double answer = 0;
    answer = (-k*pow(y,.5))/(pi*pow((.66666*y),2));
    return(answer);

}

double EulerShrinkDownOG(double k) {
    double tn = 0;
    double yn = y0og;
    int i = 1;
    while(yn > .01) {
        tn = tn +h;
        yn = yn + h*(Yprime(yn,k));        
    }
    return(tn);
    
}

double SolvableFunction(double k) {
    double answer = 0;
    double minus = 19.84;
    answer = pow(pow((EulerShrinkDownOG(k) - minus),2),.5);
    return(answer);

}

double SolvableFunctionDeriv(double k) {
    double answer = 0;
    answer = SolvableFunction(k + h)-SolvableFunction(k-h);
    answer = answer/(2*h);
    return(answer);
}

double NewtonOnSolve(double k0) {
    double kn = k0;
    while (SolvableFunction(kn)>.001) {
        kn = kn - (SolvableFunction(kn)/SolvableFunctionDeriv(kn));
        cout << "kn: " << kn << "   Time: " << SolvableFunction(kn) << endl;
    }
    cout << "kn: " <<kn << endl;
    return(kn);

}



double YPrimeN(double y, double k) {
    double answer = 0;
    answer = (-k*pow(y,.5))/(pi*pow(((-.666*y)+12),2));
    return(answer);
}

double EulerShrinkDownNew(double k) {
    double tn = 0;
    double yn = 1.985;
    int i = 1;
    while(yn > .01) {
        tn = tn +h;
        yn = yn + h*(YPrimeN(yn,k));
        
    }
    return(tn);
    
}

int main() {
    double answer = 0;
    answer = EulerShrinkDownNew(NewtonOnSolve(kNaught));
    cout << "answer: "<<answer << endl;
}

/*

This problem I solved by firstly solving for the k coeffecient using Eulers
 method to get test results for each kn and Newtons zeros method to progress though the kns
  based on an arbitrary starting volume, then I used that K and the sam
e arbitrary starting volume but if the cone turned upside down and that also required me to reverse
A(y) equation

I got the answer 77.894


*/