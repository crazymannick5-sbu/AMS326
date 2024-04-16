#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
using namespace std; 

double n1 = 10000;
double n2 = 1000000;
double pi = 3.14159265359;
double e = 2.71828182846;

double funct(double x) {
    return(sin(x)/x);
}

double rectRule(double a, double b, double n) {
    double ans = 0;
    double h = 0;
    h = (b-a)/n;
    
    for (double i = a; i<=b; i=i+h) {
        ans = ans + funct(i)*h;
    }
    return(ans);
}

double trapRule(double a, double b, double n) {
    double ans = 0;
    double h = 0;
    h = (b-a)/n;
    
    for (double i = a; (i+h)<=b; i=i+h) {
        ans = ans + (((funct(i)+funct(i+h))/2)*(h) );
    }
    return(ans);
}

double randomDouble(double lower_bound, double upper_bound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    return dist(gen);

}
double normDouble(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double Gau1Rule(double a, double b, double n) {
    double ans = 0;
    double h = 0;
    double sum = 0;
    for(int i =0; i<n; i++) {
        sum = sum+funct(randomDouble(a,b));
    }
    ans = ((b-a)/n)*sum;
    return(ans);
    
}

double Gaufunct(double x, double u, double s){
    double ans = 0;
    ans = (1.0/(s*pow((2.0*pi),.5)))*pow(e,(  (-.5)*pow(    (x-u)/s    ,2.0)    ));
    return(ans);
}

double Gau2Rule(double a, double b, double u, double s, double n) {
    double ans = 0;
    double h = 0;
    double sum = 0;
    double xi = 0;
    for(int i =0; i<n; i++) {
        xi = normDouble(u,s);
        sum = sum+((funct(xi))/Gaufunct(xi,u,s))    ;
    }
    ans = sum/n;
    return(ans);
    
}

int main() {
    cout << "rect Rule: " << rectRule(-4*pi,4*pi,n1) << endl;
    cout << "trap Rule: " << trapRule(-4*pi,4*pi,n1) << endl;
    cout << "gau1 Rule: " << Gau1Rule(-4*pi,4*pi,n2) << endl;
    cout << "gau2 Rule: " << Gau2Rule(-4*pi,4*pi,0,pi,n2) << endl;
}

/*rect Rule: 2.98432
trap Rule: 2.98432
gau1 Rule: 2.98563
gau2 Rule: 2.86402*/

//one of my results

//compile with g++ Midterm2.cpp -o a.exe
//run with a.exe