#include <iostream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <random>
using namespace std; 

double distanceFunction(double xval, double yval, double x, double y) {
    return(pow((pow((xval-x) ,2)+pow((yval-y),2)),.5));
}

double absolute(double I) {
    return(pow((pow(I,2)),.5));
}

double tophalf(double x) {
    return(sqrt(2-pow(x,2))+sqrt(absolute(x)) );
}

double bothalf(double x) {
    return(-1*(sqrt(2-pow(x,2)))+sqrt(absolute(x)) );
}
/// These are parameters to be edited/////////////////
double xmin = -1*pow((2),.5); 
double xmax = 1*pow((2),.5);
double h = .01; //Step size for heart
double hInt = .01; //Step size for area integral of heart
double hD = .01; //step size for the derivative of the heart 
double hI = .01; //step size for the integral of the perimeter arc length
//////////////////////////////////////////////////////
double radiusCheck(double xval, double yval) {
    double radmax = 0;
    double x = 0;
    double trad = 0;
    double brad = 0;
    x = xmin;

    
    while (x <= xmax) {
        double yt = tophalf(x);
        double yb = bothalf(x);
        trad = distanceFunction(xval,yval, x, yt);
        brad = distanceFunction(xval,yval,x,yb);
        if (max(brad,trad) >= radmax)  {
            radmax = max(brad,trad);
        }
        x = x+h;
    }
    return(radmax);
}


double dxdyGen(){
    double radmin = 10000;
    double stowx = 0;
    double stowy = 0;
    double checker = 0;
    double xval = xmin;
    double percent = 0;
    while(xval <= xmax) {
        cout << "x: " << xval << endl;
        /*double frac = absolute();
        if(frac >= percent + .0001 ) {
            percent = frac;
            double perc = percent * 100;
            cout << perc << "%" << endl;
        }*/
        double topval = 0;
        double botval = 0;
        topval = tophalf(xval);
        botval = bothalf(xval);
        double yval = botval;
        while(yval <= topval) {
            checker = radiusCheck(xval, yval);
            if(checker < radmin ) {
                radmin = checker;
                stowx = xval;
                stowy = yval;
            }
            yval = yval + h;
        }
        xval = xval+h;
    }
    std::cout << "Circle Radius: " << radmin << endl;
    std::cout << "x Cord Offset: " << stowx << endl;
    std::cout << "y Cord Offset: " << stowy << endl;
    return(radmin);
}


void swapRows(double** A, int row1, int row2, int N) {
    for (int i = 0; i < N; ++i) {
        double temp = A[row1][i];
        A[row1][i] = A[row2][i];
        A[row2][i] = temp;
    }
}

double tophalfINT(double x) {
    double boarder = 0;
    boarder = tophalf(xmax-.0001);
    /*if(tophalf(xmin) == tophalf(xmax)) {
        boarder = tophalf(xmax);
    } else { cout << "NOPE" << endl; exit(3); }*/
    return(tophalf(x) - boarder);

}

double bothalfINT(double x) {
    double boarder = 0;
    /*if(bothalf(xmin) == bothalf(xmax)) {
        boarder = bothalf(xmax);
    } else { cout << "NOPE" << endl;exit(3); }*/
    boarder = bothalf(xmax-.0001);
    return(bothalf(x) - boarder);
}


double simpsonsOTTOP(double a, double b, double hInt) {

    double nStart = (b - a)/h;
    int n = (int)nStart;
    double ns = (double)n;
    hInt = (b-a)/ns;

    double SUM = tophalfINT(a) + tophalfINT(b);
    for (int i = 1; i<=(n-1); i=i+2) {
        double f = i/n;
        double additive = 4*tophalfINT(a+i*hInt);
        SUM = SUM + 4*tophalfINT(a+i*hInt);
    }
    for(int i = 2; i<=(n-2); i=i+2){
        SUM = SUM + 2*tophalfINT(a+i*hInt);
    }
    return(SUM*(h/3));
}

double simpsonsOTBOT(double a, double b, double hInt) {
    double nStart = (b - a)/h;
    int n = (int)nStart;
    double ns = (double)n;
    hInt = (b-a)/ns;
    double SUM = bothalfINT(a) + bothalfINT(b);
    for (int i = 1; i<=(n-1); i=i+2) {
        double f = i/n;
        double additive = 4*bothalfINT(a+i*hInt);
        SUM = SUM + 4*bothalfINT(a+i*hInt);
    }
    for(int i = 2; i<=(n-2); i=i+2){
        SUM = SUM + 2*bothalfINT(a+i*hInt);
    }
    return(SUM*(h/3));
}

double derivativeTop(double x) {
    //std::cout << "deriv Top: " << ((tophalf(x+hD)-tophalf(x-hD))/(2*hD)) << endl;
    //std::cout << "plugged in values: " << x+hD << ", " << x-hD << endl;
    if(x+hD <= xmin) {x = x + .001;}
    if(x+hD >= xmax) {x = x -.001;}
    if(x-hD <= xmin) {x = x + .001;}
    if(x-hD >= xmax) {x = x -.001;}
    return((tophalf(x+hD)-tophalf(x-hD))/(2*hD));
}
double derivativeBot(double x) {
    //std::cout << "deriv Bot: " << ((bothalf(x+hD)-bothalf(x-hD))/(2*hD)) << endl;
    //std::cout << "plugged in values: " << x+hD << ", " << x-hD << endl;
    if(x+hD <= xmin) {x = x + .001;}
    if(x+hD >= xmax) {x = x -.001;}
    if(x-hD <= xmin) {x = x + .001;}
    if(x-hD >= xmax) {x = x -.001;}
    return((bothalf(x+hD)-bothalf(x-hD))/(2*hD));
}
double LtopIntegrand(double x) {
    return(pow((1+pow(derivativeTop(x),2)),.5));
}
double LbotIntegrand(double x) {
    return(pow((1+pow(derivativeBot(x),2)),.5));
}

double PathIntegralTop(double a, double b) {
    double N_one = (b-a)/hI;
    double N_two = (int)N_one;
    double N_three = (double)N_two;
    if ((int)N_three % 2 != 0) {
        N_three = N_three + 1;
    }
    double hDR = (b-a)/N_three;
    double SUM = LtopIntegrand(a) + LtopIntegrand(b);
    cout << "First Sum: " << SUM << endl;
    for (int i=1; i<=N_three-1; i = i+2) {
        SUM = SUM + 4*LtopIntegrand(a+i*hDR);
    }

    for (int i=2; i<=N_three-2; i =i+2) {
        SUM = SUM + 2*LtopIntegrand(a+i*hDR);
    }
    return(SUM*(hDR/3));
}

double PathIntegralBot(double a, double b) {
    double N_one = (b-a)/hI;
    double N_two = (int)N_one;
    double N_three = (double)N_two;
    if ((int)N_three % 2 != 0) {
        N_three = N_three + 1;
    }
    double hDR = (b-a)/N_three;
    double SUM = LbotIntegrand(a) + LbotIntegrand(b);
    for (int i=1; i<=N_three-1; i = i+2) {
        SUM = SUM + 4*LbotIntegrand(a+i*hDR);
    }

    for (int i=2; i<=N_three-2; i =i+2) {
        SUM = SUM + 2*LbotIntegrand(a+i*hDR);
    }
    return(SUM*(hDR/3));
}


double randomDouble(double lower_bound, double upper_bound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    return dist(gen);

}

void print2DArray(double** array, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::cout << array[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void print1DArray(double* array, int size) {
    for (int i = 0; i < size; ++i) {
        std::cout << array[i] << " ";
    }
    std::cout << std::endl;
}

double** ArrayInstant(double up, double down, int N) {
    double **Main = new double*[N];
    for (int k = 0; k<N; k++) {
        Main[k] = new double[N];
    }
    for (int i = 0; i < N; ++i) {
        for(int j = 0; j<N; ++j) {
            Main[i][j] = randomDouble(up, down);
        }
    }
    print2DArray(Main, N);
    return(Main);

}


double* MatVectMult(double **A, double *X, int N) {
    double* result = new double[N];

    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for(int j = 0; j < N; ++j) {
            result[i] = result[i] + (A[i][j]*X[j]);
        }
    } 
    return(result);
}




double** makeDiagDominant(double** A, double* b, int N) {
    for (int i = 0; i < N; ++i) {
        double diagonal = A[i][i];
        double sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (j != i) {
                sum += abs(A[i][j]);
            }
        }
        if (abs(diagonal) <= sum) {
            int maxIndex = i;
            double maxDiagonal = abs(diagonal);
            for (int k = i + 1; k < N; ++k) {
                double currDiagonal = abs(A[k][k]);
                if (currDiagonal > maxDiagonal) {
                    maxDiagonal = currDiagonal;
                    maxIndex = k;
                }
            }
            swapRows(A, i, maxIndex, N);
            double temp = b[i];
            b[i] = b[maxIndex];
            b[maxIndex] = temp;
        }
    }
    double** V = new double*[N];
    for(int i = 0; i<N; i++) {
        V[i] = new double[N+1];
    }
    for(int i = 0; i<N; i++) {
        for(int j = 0; j<N; j++) {
            V[i][j] = A[i][j];
        }
        V[i][N] = b[i];
    }
    return(V);
}







double* Solver(double **A, int N, int n) {
    cout << "check 3" << endl;
    double** D = new double*[N];
    double** D_1 = new double*[N];
    double** R = new double*[N];
    double* b = new double[N];
    double* X = new double[N];
    for(int i = 0; i<N; i++) {
        D[i] = new double[N];
        D_1[i] = new double[N];
        R[i] = new double[N];
    }
    cout << "check 2" << endl;
    for(int i = 0; i < N; i++) {
        b[i] = 1;
        X[i] = 0;
        for(int j = 0; j < N; j++) {
            D[i][j] = 0;
            D_1[i][j] = 0;
            R[i][j] = A[i][j];
        }
    }
    cout << "check 1" << endl;
    for(int i = 0; i < N; i++) {
        D[i][i] = A[i][i];
        D_1[i][i] = 1/(D[i][i]);
        R[i][i] = 0; 
    }
    cout << "X Array: " << endl;
    print1DArray(X,N);
    cout << "b Array: " << endl;
    print1DArray(b,N);
    cout << "D Array: " << endl;
    print2DArray(D,N);
    cout << "D_1 Array: " << endl;
    print2DArray(D_1,N);
    cout << "R Array: " << endl;
    print2DArray(R,N);
    for(int k=0; k<=n; k++) {
        double* U = MatVectMult(R,X,N);
        double* L = new double[N];
        for(int ii=0; ii<N; ++ii) {
            L[ii] = b[ii]-U[ii];
        }
        X = MatVectMult(D_1,L,N);
    }
    print1DArray(X, N);
    return(X);
}


int main() {

    double radmin = dxdyGen();
    double AOC = 3.1415926 * pow(radmin,2);
    std::cout << "Area of Circle: " << AOC << endl;
    double AOH = absolute(simpsonsOTTOP((xmin+.0001),(xmax-.0001),hInt) - simpsonsOTBOT((xmin+.0001),(xmax-.0001),hInt));
    std::cout << "Area of Heart: " << AOH << endl;
    std::cout << "Area left Around the Heart : " << AOC-AOH << endl;
    double COH = absolute(PathIntegralTop((xmin+.01),(xmax-.01))) + absolute(PathIntegralBot((xmin+.01),(xmax-.01)));
    std::cout << "Circumference of the Heart: " << COH << endl;
    int N = 6;
    std::cout << "check 7" << endl;
    double *b = new double[N];
    for(int i=0; i<N ; i++) {b[N] = 1;}
    double **A = ArrayInstant(15.0,1.0,N);
    double **U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    U = makeDiagDominant(A,b,N);
    for(int i =0; i<N; i++) {
        for(int j=0; j<N; j++) {
            A[i][j]=U[i][j];
        }
    }
    std::cout << "DIAGONALIZED MATRIX: " << endl;
    print2DArray(A,N);
    std::cout << "check 4" << endl;
    double *X = Solver(A, N, 30); 
    std::cout << "check 5" << endl;
    print1DArray(MatVectMult(A,X,N), N);
    std::cout << "check 6" << endl;
    std::cout<<"FINISHED!" << endl;
    
    
    
    
    
    /*double** T = new double*[4];
    double* T_1 = new double[4];
    for(int i = 0; i<N; i++) {
        T[i] = new double[4];
        T_1[i] = 2;
        for (int j = 0; j<4; j++) {
            T[i][j] = j+i;
        }
    }
    std::cout << "T Matrix: " << endl;
    print2DArray(T,4);
    std::cout << "T_1 Vector: " << endl;
    print1DArray(T_1,4);
    std::cout << "solve: " << endl;
    print1DArray(MatVectMult(T,T_1,4),4);*/

    
    
}