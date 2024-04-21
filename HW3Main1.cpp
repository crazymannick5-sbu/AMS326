#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
using namespace std; 

double h = .0001;
double xBig = .7693;
double yBig = .7693;
double zeroError = .00001;
int number = 0;
double pi = 3.14159265359;
int sampleSize = 10000;

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

double randomDoubleBetween(double min, double max) {
    std::random_device rd; // Obtain a random seed from the OS entropy device
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(min, max); // Define the range

    return dis(gen);
}

double fLxy(double x, double y) {
    double answer = 0;
    answer = pow((pow(x,2)+pow(y,2)),3)-(4*pow(x,2)*pow(y,2));
    return(answer);
}
double fLdx(double x, double y) {
    double answer = 0;
    answer = (fLxy((x+h),(y))-fLxy((x-h),(y)))/(2*h);
    return(answer);
}
double fLdy(double y, double x) {
    double answer = 0;
    answer = (fLxy((x),(y+h))-fLxy((x),(y-h)))/(2*h);
    return(answer);
}
double *fLx(double y) {
    double xCheck = -1 * xBig;
    int numb = 0;
    double check  = fLxy(xCheck,y)/(pow(pow(fLxy(xCheck,y),2),.5));

    while(fLxy(xCheck,y) < xBig) {
        if(fLxy(xCheck,y) <0 && check > 0) {
            check = -1;
            numb = numb + 1;
        }
        if(fLxy(xCheck,y) >0 && check <0) {
            check = 1;
            numb = numb + 1; 
        }
        xCheck = xCheck +h; 
    }
    number = numb; 
    double* answer = new double[numb];
    double* initial = new double[numb];

    xCheck = -1 * xBig;
    numb = 0;
    check  = fLxy(xCheck,y)/(pow(pow(fLxy(xCheck,y),2),.5));
    while(fLxy(xCheck,y) < xBig) {
        if(fLxy(xCheck,y) <0 && check > 0) {
            check = -1;
            initial[numb] = xCheck;
            numb = numb + 1;
        }
        if(fLxy(xCheck,y) >0 && check <0) {
            check = 1;
            initial[numb] = xCheck;
            numb = numb + 1; 
        }
        xCheck = xCheck +h; 
    }
    

    double e = 1;
    double xi = 1;
    int count = 0;
    while(count < numb) {
        e = 1;
        xi = initial[count];
        //cout << "checking on " << count << ": " << initial[count] << endl;

        while(e > zeroError) {
            xi = xi - fLxy(xi,y)/fLdx(xi,y);
            e = pow(pow((  0-fLxy(xi,y)  ),2),.5);

        }
        answer[count] = xi; 
        count = count + 1;
    }

    
    return(answer);
}
double *fLy(double x) {
    double yCheck = -1 * yBig;
    int numb = 0;
    double check  = fLxy(yCheck,x)/(pow(pow(fLxy(yCheck,x),2),.5));

    while(fLxy(yCheck,x) < yBig) {
        if(fLxy(yCheck,x) <0 && check > 0) {
            check = -1;
            numb = numb + 1;
        }
        if(fLxy(yCheck,x) >0 && check <0) {
            check = 1;
            numb = numb + 1; 
        }
        yCheck = yCheck +h; 
    }
    number = numb; 
    double* answer = new double[numb];
    double* initial = new double[numb];

    yCheck = -1 * yBig;
    numb = 0;
    check  = fLxy(yCheck,x)/(pow(pow(fLxy(yCheck,x),2),.5));
    while(fLxy(yCheck,x) < yBig) {
        if(fLxy(yCheck,x) <0 && check > 0) {
            check = -1;
            initial[numb] = yCheck;
            numb = numb + 1;
        }
        if(fLxy(yCheck,x) >0 && check <0) {
            check = 1;
            initial[numb] = yCheck;
            numb = numb + 1; 
        }
        yCheck = yCheck +h; 
    }
    

    double e = 1;
    double yi = 1;
    int count = 0;
    while(count < numb) {
        e = 1;
        yi = initial[count];
        //cout << "checking on " << count << ": " << initial[count] << endl;

        while(e > zeroError) {
            yi = yi - fLxy(yi,x)/fLdx(yi,x);
            e = pow(pow((  0-fLxy(yi,x)  ),2),.5);

        }
        answer[count] = yi; 
        count = count + 1;
    }

    
    return(answer);
}

int cloverCheck(double x, double y) {
    bool answer = 0;
    double* outputs = fLy(x);
    int moreThanCount = 0;
    int len =number;
    //cout << "length " << len << endl;
    for (int i = 0; i< len; i++) {
        if (y > outputs[i]) {
            moreThanCount = moreThanCount + 1;
            //cout << "y is more than " << outputs[i] << endl;
        }
    }
    if(moreThanCount%2 != 0) {
        answer = 1;

    }
    return(answer);

}

int ReCheck(double y1, double x1, double H, double x, double y ) {
    int ans = 0;
    double ya = tan(H)*(x-x1) + (.5/cos(H)) + (y1);
    double yb = tan(H)*(x-x1) - (.5/cos(H)) + (y1);
    double yc = tan(H+(pi/2))*(x-x1) - (((1/pow(2,.5))/2)/cos(H+(pi/2))) + (y1);
    double yd = tan(H+(pi/2))*(x-x1) + (((1/pow(2,.5))/2)/cos(H+(pi/2))) + (y1);


    /*cout << "a" << ya << endl;
    cout << "b" << yb << endl;
    cout << "c" << yc << endl;
    cout << "d" << yd << endl;*/
    int check1 = 0;
    int check2 = 0;

    if((y>ya&&y<yb) || (y<ya&&y>yb)) {
        check1 = 1;
        //cout << "one passed" << endl;
    }
    if((y>yc&&y<yd) || (y<yc&&y>yd)) {
        //cout << "two passed" << endl;
        check2 = 1;
    }

    if(check2 == 1 && check1 == 1) {
        
        ans = 1;
    }

    return(ans);
}

int BothCheck(double x, double y, double y1, double x1, double H) {
    int ans = 0;
    if (cloverCheck(x,y)==1 && ReCheck(y1,x1,H,x,y)==1) {
        ans =1;
    }
    return(ans);
}

double AreaFunction(double x1, double y1, double H) {
    double i = 0;
    double area = 0;
    double areaWhole = 4*xBig*yBig; 
    double areaRatio = 0;
    double yesSum= 0;
    for (i=0; i<sampleSize; i++) {
        double x = randomDoubleBetween(-1*xBig,xBig);
        double y = randomDoubleBetween(-1*yBig,yBig);
        yesSum = yesSum + BothCheck(x,y,y1,x1,H);
        showProgressBar(i,sampleSize);
    }
    areaRatio = yesSum/i;
    area = areaWhole * areaRatio;
    return(area);
}


int main() {
    double * answer = fLx(0.3);
    cout << "1: " << answer[0] << endl;
    cout << "2: " << answer[1] << endl;
    cout << "3: " << answer[2] << endl;
    cout << "4: " << answer[3] << endl;
    double * answer1 = fLy(0.4321);
    cout << "1: " << answer1[0] << endl;
    cout << "2: " << answer1[1] << endl;
    cout << "3: " << answer1[2] << endl;
    cout << "4: " << answer1[3] << endl;
    int ans = cloverCheck(-1, 1);
    cout << "Checking: " << ans << endl;
    int ans1 = ReCheck(0,-.2,-.54,-.5,0);
    cout << "Rect Checking: " << ans1 << endl;
    double area = AreaFunction(-.9,.7,.5);
    cout << "area calculation: " << area << endl;
}