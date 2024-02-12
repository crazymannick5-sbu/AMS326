#include <iostream>
#include <cmath>
#include <functional>


/*double fFunction(double x){
    double A = 0;
    double B = 0;
    B=-1*pow(x,2.0);
    A = exp(B);
    return(A-pow(x,3.0));

}*/


std::function<double(double)> calculateDerivative(const std::function<double(double)>& func, double epsilon = 1e-6) {
    return [func, epsilon](double x) {
        return (func(x + epsilon) - func(x - epsilon)) / (2 * epsilon);
    };
}

double fFunction(double x) {
    return(pow(x,3.0)+x-1);
}

double BisectionMain(double a1, double b1) {
    double Error = .00001;
    /*std::cout << "a1: " << a1 << std::endl;*/
    /*std::cout << "b1: " << b1 << std::endl;*/
    double x1 = (a1 + b1)/2.0;
    /*std::cout << "x1: " << x1 << std::endl;*/
    double fa1 = fFunction(a1);
    double fx1 = fFunction(x1);
    double Mult = fa1 * fx1;
    /*std::cout << "M: " << Mult << std::endl;*/
    double ans = 0;
    
    if(pow(pow(a1-b1,2),.5) < Error) {ans = (a1+b1)/2; std::cout << ans << std::endl;    } else {
        if(Mult < 0) {
            ans = BisectionMain(a1, x1);
        } else {
            ans = BisectionMain(x1,b1);
        }
    }
    return(ans);
    /*I chose to perform this part of the assignment in a recursive manner instead of with a for loop,
    there are of course ups and downs to each method but this wasn't a choice made on the the effeciency
    of the program and rather my interest in the code. */
}

double NewtonMain(double x1){
    /*This is my newtonian method function. This function relies on the derivative of a function. I got a
    little side tracked and decided to research some numerical methods of derivation and wrote a function
    that calculates the derivative instead of relying on user based code changes. */
    auto dF = calculateDerivative(fFunction);
    double e = .0000001;
    double maxI = 5000;
    double ans = 0;
    for(int i = 0; i<maxI; i++) {
        x1 = x1 - fFunction(x1)/dF(x1);
        if(pow(pow(fFunction(x1),2),.5) < e ) {
            ans = x1;
            break;
        }
    }
    std::cout << x1 << std::endl;
    return(ans);


}

double SCantMain(double x1, double x0){
    /* */
    double x2 = 0.0;

    double e = .0000001;
    double maxI = 5000;
    double ans = 0;
    for(int i = 0; i<maxI; i++) {
        x2 = x1 - ((x1-x0)*((fFunction(x1))/(fFunction(x1)-fFunction(x0))));
        x0 = x1;
        x1 = x2;
        if(pow(pow(fFunction(x1),2),.5) < e ) {
            ans = x1;
            break;
        }
    }
    std::cout << x1 << std::endl;
    return(ans);


}



int main() {
    double answer = 0;
    answer = BisectionMain(0,1);
    double answer1 = 0;
    answer1 = NewtonMain(.4);
    double answer3 = 0;
    answer3 = SCantMain(0,1);
 
}