#ifndef __MAIN_H__
#define __MAIN_H__
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry> //to use quaternions
#include <string.h>
#include <cmath>
using namespace Eigen;
using namespace std;
// In this file you can find all the useful and common functions and constants used

// vectors of the A distance (metres), D distances (metres) and angle variables
const VectorXd A = (VectorXd(6) << 0, -0.425, -0.3922, 0, 0, 0).finished();
const VectorXd D = (VectorXd(6) << 0.1625, 0, 0, 0.1333, 0.0997, 0.0996).finished();
const VectorXd ALPHA = (VectorXd(6) << M_PI / 2, 0, 0, M_PI / 2, -M_PI / 2, 0).finished();

// Compute transformation matrix from i to j joint
Matrix4d computeTransformationMatrix(double th, double alpha, double distance, double a)
{
    Matrix4d Tij;
    Tij << cos(th), -sin(th) * cos(alpha), sin(th) * sin(alpha), a * cos(th),
        sin(th), cos(th) * cos(alpha), -cos(th) * sin(alpha), a * sin(th),
        0, sin(alpha), cos(alpha), distance,
        0, 0, 0, 1;
    return Tij;
}

// evaluates if a number is very close to zero
bool almostZero(double x)
{
    return (abs(x) < 1e-7);
}

int print_eigen(string str, Matrix4d m)
{
    cout << str << endl
         << m << endl;
    return 0;
}

#endif