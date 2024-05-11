#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
#include <complex>
using namespace std; 
double e = 2.71828182;
double pi = 3.141592659;
double XMax = 4;
double YMax = 1;
double XMin = -4;
double YMin = 0;
int N = 10000;
std::complex<double> i(0,1);


double randomDoubleBetween(double min, double max) {
    std::random_device rd; // Obtain a random seed from the OS entropy device
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(min, max); // Define the range

    return dis(gen);
}

double SBlue(double x) {
    double answer = 0;
    answer = 1/(1+pow(e,(-1*x)));
    return(answer);
}

double SGold(double x) {
    double answer = 0;
    answer = 1/(1+pow(e,(-1*x)));
    return(answer);
}

bool CheckGold(double xi, double yi, double H, double d) {
    bool answer = false;
    double xf = xi + d*cos(H);
    double yf = yi + d*sin(H);

    int a1 = 0;
    int a2 = 0;

    if(SGold(xi) < yi && SGold(xf) > yf) {
        a1 = 1;
    }

    if(SGold(xi) > yi && SGold(xf) < yf) {
        a2 = 1;
    }

    if(a1 ==1 || a2 == 1) {
        answer = true;
    }

    return(answer);
    

}
bool CheckBlue(double xi, double yi, double H, double d) {
    int answer = false;
    double xf = xi + d*cos(H);
    double yf = yi + d*sin(H);

    int a1 = 0;
    int a2 = 0;

    if(SBlue(xi) < yi && SBlue(xf) > yf) {
        a1 = 1;
    }

    if(SBlue(xi) > yi && SBlue(xf) < yf) {
        a2 = 1;
    }

    if(a1 ==1 || a2 == 1) {
        answer = true;
    }
    
    return(answer);
}

int Check(double xi, double yi, double H, double d) {


    int answer = 0;
    if(CheckGold(xi,yi,H,d) || CheckBlue(xi,yi,H,d)) {
        answer = 1;
    }
    return(answer);
}

bool CheckPoint(double xi,double yi,double H,double d) {
    bool answer = true;
    double xf = xi + d*cos(H);
    double yf = yi + d*sin(H);
    //cout << "xf: " << xf << endl;
    ////cout << "yf: " << yf << endl;
    if(xf > 4 || xf < -4 || yf > 1 || yf < 0) {
        answer = false;
        ///cout << "fail try again " << endl;
        
    } 

    return(answer);
}

void showProgressBar(int progress, int total) {
    const int barWidth = 50;

    float percentage = (float)progress / total;
    int numBars = barWidth * percentage;

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < numBars)
            std::cout << "=";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(percentage * 100.0) << "%\r";
    std::cout.flush();
}

void NumericalCheck(double d) {
    int sum = 0;
    double xi = randomDoubleBetween(-4,4);
    double yi = randomDoubleBetween(0,1);
    double H = randomDoubleBetween(0,(2*pi));
    bool Checker = false;
    //cout << "start" << endl;
    for(int i = 0; i < N; i++) {
        ////cout << "for loop: " << i << endl;
        double xi = randomDoubleBetween(-4,4);
        double yi = randomDoubleBetween(0,1);
        double H = randomDoubleBetween(0,(2*pi));
        bool Checker = false;
        //cout << "for loop" << endl;
        while (!Checker) {

            //cout << "while" << endl;
            //cout <<"xi: " << xi << endl;
            //cout <<"yi: " << yi << endl;
            //cout << "H: " << H <<endl;
            if (!CheckPoint(xi,yi,H,d)) {
                xi = randomDoubleBetween(-4,4);
                
                yi = randomDoubleBetween(0,1);
                H = randomDoubleBetween(0,(2*pi));
                Checker =false;
                

            } else {
                Checker = true;
                //cout << "good";
            }
        }
        sum = sum + Check(xi,yi,H,d);

        
        showProgressBar(i,N);
    }

    double ND = N;
    double sumD = sum;
    cout << (sumD*100)/ND << " percent of the needles made it in the box with d = " << d << endl;
}

int main() {
    //cout << "start" << endl;
    NumericalCheck(.25);
    NumericalCheck(.5);
    NumericalCheck(1);
    NumericalCheck(2);
}

/*The only tricky part abou this problem was deciding what to do about keeping the needles
in the box, but I ended up deciding to simply reroll the needle if it was going to land 
outside the box. */