#ifndef __MOTION_PLANNER_H__
#define __MOTION_PLANNER_H__
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry> //to use quaternions
#include <string.h>
#include <cmath>
// #include "ros/ros.h"
using namespace Eigen;
using namespace std;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 8, 1> Vector8d;
typedef Matrix<double, 3, 1> Vector3d;
typedef Matrix<double, Dynamic, 8> Path;

// vectors of the A distance (metres), D distances (metres) and angle variables
const VectorXd A = (VectorXd(6) << 0, -0.425, -0.3922, 0, 0, 0).finished();
const VectorXd D = (VectorXd(6) << 0.1625, 0, 0, 0.1333, 0.0997, 0.0996).finished();
const VectorXd ALPHA = (VectorXd(6) << M_PI / 2, 0, 0, M_PI / 2, -M_PI / 2, 0).finished();
const VectorXd a = (VectorXd(6) << 0, -0.425, -0.3922, 0, 0, 0).finished();
const VectorXd d = (VectorXd(6) << 0.1625, 0, 0, 0.1333, 0.0997, 0.0996).finished();
const VectorXd alpha = (VectorXd(6) << M_PI / 2, 0, 0, M_PI / 2, -M_PI / 2, 0).finished();
// useful deltas: used in interpolation-->time and path duration
const double dt = 0.1;
const double path_dt = 10.0;

Matrix6d jacobian(Vector6d js);
// Compute transformation matrix from i to j joint
Matrix4d computeTransformationMatrix(double th, double alpha, double distance, double a);

// evaluates if a number is very close to zero
bool almostZero(double x);
Matrix4d rotation_180z_axis_and_offset();
int print_eigen(string str, MatrixXd m);

// Functio to compute the jacobian matrix
Matrix6d computeJacobian(Vector6d th);

// the function directKin returns the direct Kinematic esxpressed in a matrix Tm (Transformation Matrix)
Matrix4d directKin(VectorXd Th);

// function to compute the inverse kinematics:known the final pose we want to know the joint angles
Matrix<double, 6, 8> inverseKin(Vector3d p60, Matrix3d R60, double scaleFactor);

// to interpolate two vectors and correctly compute the SLERP,
// first of all we need to interpolate the 2 corresponding points in the space.
// We find an "interpolated" point, which we can call it p_i (Point of Interpolation),
// since it is in the space, it will be a 3D vector.
Vector3d lerp(double time, Vector3d p1, Vector3d p2);

// as we did for position p1 and p2, we need to do the same with the respective quaternions, and find the "interpolated" quaternion
Quaterniond myslerp(double time, Quaterniond q1, Quaterniond q2);

#endif