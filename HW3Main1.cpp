#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
#include <fstream>
#include <string>



using namespace std; 

double h = .0001;
double xBig = .7693;
double yBig = .7693;
double zeroError = .00001;
int number = 0;
double pi = 3.14159265359;
int sampleSize = 1000;
int sampleSize2= 100000;
double w = 25;
double v0 = 100;
double k = w/v0;

void runPythonScript(const std::string& pythonFileName) {
    std::string command = "python " + pythonFileName;
    int result = std::system(command.c_str());
    if (result == 0) {
        std::cout << "Python script executed successfully." << std::endl;
    } else {
        std::cerr << "Error: Failed to execute Python script." << std::endl;
    }
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
void showProgressBar2(double progress, double total) {
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
        //showProgressBar(i,sampleSize);
    }
    areaRatio = yesSum/i;
    area = areaWhole * areaRatio;
    return(area);
}
double HdAreaFunction(double x1, double y1, double Hi) {
    double answer = 0;
    answer = (AreaFunction(x1,y1,(Hi+h))-AreaFunction(x1,y1,(Hi-h)))/(2*h);
    return(answer);
}
double XdAreaFunction(double xi, double y1, double H) {
    double answer = 0;
    answer = (AreaFunction((xi+h),y1,H)-AreaFunction((xi-h),y1,H))/(2*h);
    return(answer);
}
double YdAreaFunction(double x1, double yi, double H) {
    double answer = 0;
    answer = (AreaFunction(x1,(yi+h),H)-AreaFunction(x1,(yi-h),H))/(2*h);
    return(answer);
}
void gradientReturn() {
    double xn = -1;
    double yn = 1;
    double Hn = 0;
    double output0 = AreaFunction(xn,yn,Hn);
    double outputn = 0;
    double currentStepFactor = 1.0;
    double max = 1/(pow(2,.5));
    double factor = 0;
    
    double Hdg = HdAreaFunction(xn,yn,Hn);
    double Xdg = XdAreaFunction(xn,yn,Hn);
    double Ydg = YdAreaFunction(xn,yn,Hn);
    double Hdgn = 0;
    double Xdgn = 0;
    double Ydgn = 0;
    
    cout << "start" << endl;
    int dec = 1;
    while(outputn <= .75*max) {
        outputn = AreaFunction(xn,yn,Hn);
        Hdgn = HdAreaFunction(xn,yn,Hn);
        Xdgn = XdAreaFunction(xn,yn,Hn);
        Ydgn = YdAreaFunction(xn,yn,Hn);
        factor = (max - outputn)/max;
        xn = xn + factor*Xdgn;
        yn = yn + factor*Ydgn;
        Hn = Hn + factor*Hdgn;

        //showProgressBar2(outputn, .707106);
        cout << "area: " << outputn <<  endl;  
        cout << "   Hdgn: " << Hdgn << endl;
        cout << "   Xdgn: " << Xdgn << endl;
        cout << "   Ydgn: " << Ydgn << endl;

        
        
        
    }
    cout << "gradientReturn: "  <<outputn << endl;
}

double CentDet(double L, double w, double d, double cent) {
    double center = randomDoubleBetween(0,w);
    double sum = -1*center+w;
    int answer = 0;
    while(sum<= d) {
        answer = answer + 1; 
        sum = sum+w;
    }
    return(answer);
}

void printVectorOfVectors(const std::vector<std::vector<double>>& vec) {
    // Determine the maximum size of the inner vectors
    size_t maxInnerSize = 9;
    for (const auto& innerVec : vec) {
        if (innerVec.size() > maxInnerSize) {
            maxInnerSize = innerVec.size();
        }
    }

    // Print header row with column indices
    std::cout << "  ";
    for (size_t i = 0; i < maxInnerSize; ++i) {
        std::cout << i << "\t";
    }
    std::cout << std::endl;

    // Print each inner vector in tabulated form
    for (size_t row = 0; row < vec.size(); ++row) {
        std::cout << row << " ";
        for (size_t col = 0; col < vec[row].size(); ++col) {
            std::cout << vec[row][col] << "  "<< "\t";
        }
        // Pad with empty cells if necessary
        for (size_t i = vec[row].size(); i < maxInnerSize; ++i) {
            std::cout << "\t";
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<double>> ProbDetermine(double L, double w, double d) {
    double center = randomDoubleBetween(0,w);
    double max = 0;
    std::vector<double> ref;
    std::vector<double> count;
    std::vector<std::vector<double>> answer;
    
    if(d/w<=1.0) {
        max = 1;
    } else {
        if(static_cast<double>(static_cast<int>(d/w))-(d/w) != 0 ) {
            max = 1 + static_cast<int>(d/w);
        } else {
            max = static_cast<int>(d/w);
        }
    }
    for(int i = 0; i<=max; i++) {
        count.push_back(0);
        ref.push_back(i+1);
    }
    for(int i = 0; i<sampleSize2; i++) {
        center = randomDoubleBetween(0,w);
        double aN = CentDet(L,w,d,center);
        count[std::round(aN)] = count[std::round(aN)] + 1; 
    }
    cout << "max " << max << endl;
    for(int i = 0; i<=max; i++) {
        
        count[i] = count[i] / sampleSize2;
        cout<< count[i] << endl;
    }
    answer.push_back(ref);
    answer.push_back(count);
    return(answer);     
     
}

/*void writeToCSV(double num1, double num2, const string& fileName) {
    ofstream outputFile(fileName);
    if (!outputFile.is_open()) {
        cout << "Error: Unable to open file " << fileName << endl;
        return;
    }
    
    outputFile << num1 << "," << num2 << endl;
    cout << num1 << "," << num2 << endl;
    
    outputFile.close();
    //cout << "Data written to " << fileName << endl;
}*/

void deleteFile(const std::string& filename) {
    if (std::remove(filename.c_str()) != 0) {
        std::cerr << "Error deleting file: " << filename << std::endl;
        //return false;
    } else {
        std::cout << "File deleted successfully: " << filename << std::endl;
        //return true;
    }
}

void writeToCSV(double num1, double num2, const std::string& filename) {
    std::ofstream file(filename, std::ios_base::app); // Open file in append mode
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    file << num1 << "," << num2 << "\n"; // Write numbers in CSV format
    file.close(); // Close the file
}

long double yPrime(long double yn, long double tn) {
    long double e = 1;
    long double pn = 0;
    pn = (yn/tn)-k*pow((1+pow((yn/tn),2)),.5);

    return(pn);


}

/*void plotCSV(const string& fileName) {
    vector<double> xValues;
    vector<double> yValues;

    ifstream inputFile(fileName);
    if (!inputFile.is_open()) {
        cout << "Error: Unable to open file " << fileName << endl;
        return;
    }

    string line;
    while (getline(inputFile, line)) {
        stringstream ss(line);
        string token;
        getline(ss, token, ',');
        double x = stod(token);
        getline(ss, token, ',');
        double y = stod(token);
        xValues.push_back(x);
        yValues.push_back(y);
    }

    inputFile.close();

    plt::plot(xValues, yValues);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::title("Data Plot");
    plt::show();
}*/


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

void ReverseEulerEstimate(long double y0, long double x0, long double hE, long double n) {
    long double yn = y0;
    long double xn = 0;
    for(long double i = 1.0; i <= n; i = i+1) {
        xn = x0-(hE*i);
        yn = yn - (hE*yPrime(yn,xn));
        long double xIN = log2(xn);
        if(trunc(xIN)==xIN) {
            //cout << "log2(tn) = " << xIN << "   y: " << yn << endl;
        }
        writeToCSV(xn,yn,"Data.csv");
        showProgressBar(i,n);

        
    }
    runPythonScript("Plot.py");
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
    /*double * answer = fLx(0.3);
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
    //double area = AreaFunction(-1,1,-1.570696);
    //cout << "area calculation: " << area << endl;
    gradientReturn();*/

    //FIRST ASSIGNMENT HERE;

    std::vector<std::vector<double>> Output; 
    Output.push_back(ProbDetermine(1,1,.1)[1]);
    Output.push_back(ProbDetermine(1,1,.5)[1]);
    Output.push_back(ProbDetermine(1,1,1)[1]);
    Output.push_back(ProbDetermine(1,1,2)[1]);
    Output.push_back(ProbDetermine(1,1,3)[1]);
    printVectorOfVectors(Output);


    deleteFile("Data.csv");
    ReverseEulerEstimate(0,160,.01,15999);
    
    
}