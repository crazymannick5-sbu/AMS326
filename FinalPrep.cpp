#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
#include <complex>
using namespace std; 

double e = 2.71828182;
double pi = 3.141592659;
std::complex<double> i(0,1);

std::complex<double> pruW(int n) {
    std::complex<double> exponent;
    std::complex<double> answer;
    exponent = -1.0*i*(2*(pi/n));
    return(pow(e,exponent));

}

std::complex<double> yKElement(int n, int k, const std::vector<double>& X) {
    std::complex<double> answer=0;
    std::complex<double> sum = 0;

    for(int i = 0; i <n; i++) {
        sum = sum + (X[i]*pow(pruW(n),(i*k)));

    }
    answer = pow(n,-.5)*sum;
    return(answer);
    
}

std::vector<std::complex<double>> YVectorForm(int n, const std::vector<double> & X) {
    std::vector<std::complex<double>> answer;
    for(int i = 0; i<n; i++) {
        answer.push_back(yKElement(n,i,X));
    }
    return(answer);
}


