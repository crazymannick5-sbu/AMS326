#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

struct Point {
    double x;
    double y;
};

std::vector<double> lagrangeInterpolation(const std::vector<Point>& points) {
    std::vector<double> coefficients(points.size(), 0.0);

    for (size_t i = 0; i < points.size(); ++i) {
        double term_coefficient = points[i].y;
        for (size_t j = 0; j < points.size(); ++j) {
            if (i != j) {
                term_coefficient /= (points[i].x - points[j].x);
                coefficients[i] += term_coefficient;
                coefficients[j] *= -points[i].x / (points[i].x - points[j].x);
            }
        }
    }

    return coefficients;
}

void printCoefficients(const std::vector<double>& coefficients) {
    std::cout << "Coefficients of the interpolated polynomial:" << std::endl;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        std::cout << "Coefficient for x^" << i << ": ";
        if (i == 0) {
            std::cout << coefficients[i];
        } else {
            std::cout << std::fixed << std::setprecision(2) << coefficients[i];
        }
        std::cout << std::endl;
    }
}

int main() {
    std::vector<Point> points = {{0, 1}, {2, 2}, {3, 4}}; // Test points

    std::vector<double> coefficients = lagrangeInterpolation(points);
    printCoefficients(coefficients);

    return 0;
}
