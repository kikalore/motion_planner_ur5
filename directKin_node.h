#ifndef __DIRECTKIN_NODE_H__
#define __DIRECTKIN_NODE_H__
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <string.h>
#include <cmath>
#include "main.h"
using namespace Eigen;
using namespace std;

// the function directKin returns the directKinematic esxpressed in a matrix Tm (Transformation Matrix)
void directKin(VectorXd &Th, double scaleFactor, VectorXd &pe, Matrix3d &Re, Matrix4d &Tm)
{
    // compute single transformations
    Matrix4d T10 = computeTransformationMatrix(Th(0), ALPHA(0), D(0), A(0));
    Matrix4d T21 = computeTransformationMatrix(Th(1), ALPHA(1), D(1), A(1));
    Matrix4d T32 = computeTransformationMatrix(Th(2), ALPHA(2), D(2), A(2));
    Matrix4d T43 = computeTransformationMatrix(Th(3), ALPHA(3), D(3), A(3));
    Matrix4d T54 = computeTransformationMatrix(Th(4), ALPHA(4), D(4), A(4));
    Matrix4d T65 = computeTransformationMatrix(Th(5), ALPHA(5), D(5), A(5));
    // compute transformation matrix from 0 to 6(T60=Tm)
    Tm = (T10) * (T21) * (T32) * (T43) * (T54) * (T65);
    // pe is the end effector postion: it's coincides with the first 3 rows of the last column( <3,1> vector) of Tm
    pe = Tm.block<3, 1>(0, 3); // extract from row 0 and column 3 (so 4th column) a vector of 3x1 elements
    // Re is the rotation matrix-->3x3--> first 3 rows and first 3 columns of Tm
    Re = Tm.block<3, 3>(0, 0);
}
#endif
