#include "motion_planner.h"
using namespace Eigen;
using namespace std;
//testing
int main()
{
    /*DIRECT KINEMATICS COMPUTATION OF UR5
    We have defined the followed:
    ->Th: the ur5 six joint angles
    ->pe: cartesian position of the end effector
    ->Re: rotation matrix of the end effector
    */
    // VectorXd Th(6), pe(3);
    // Matrix3d Re;
    // Matrix4d Tm;
    // Th << 0, 0, 0, 0, 0, 0;
    // double scaleFactor = 10;
    // // direct kinematics
    // directKin(Th, scaleFactor, pe, Re, Tm);
    // // Print direct kinematics matrix
    // print_eigen("DirectKinematics", Tm);
    // cout << endl;



//     /*INVERSE KINEMATICS COMPUTATION OF UR5
//     We have defined the followed:
//     ->Th: the ur5 six joint angles to FIND
//     ->pe: desidered cartesian position of the end effector
//     ->Re: desidered rotation matrix of the end effector
//     */
//    // position and orientation desidered of the end-effector
//     VectorXd pe60(3);
//     pe60 << 1, 1, 1; // to modify with the desidered position-->vision node
//     Matrix3d Re;
//     Re << 1, 1, 1,
//         1, 1, 1,
//         1, 1, 1; // to modify with the desidered orientation-->vision node

//     // joint angles to find
//     MatrixXd JointAngles(6,8);

//     // final transformation matrix
//     Matrix4d Tm;

//     // compute inverse kinematics
//     inverseKin(JointAngles,1.0, pe60, Re, Tm);

//     cout << "Possible configurations of joint angles:\n"<< JointAngles << endl; 

    //JACOBIAN
    Matrix6d myJac;
    Matrix4d transMatrix;
    Vector6d myJointVariables(6);
    myJointVariables<<-0.0001,  -0.00016, -0.00006, -0.00014, -0.00004,  0.00001;
    myJac=computeJacobian(myJointVariables);
    //transMatrix = directKin(myJointVariables);
    cout<<myJac<<"\n";

}

// Matrix6d computeJacobian(Vector6d th)
// {
//     Matrix6d jacobianMatrix;
//     Vector6d J1(6, 1); // first column of the jacobian matrix
//     J1 << D(4) * (cos(th(0)) * cos(th(4)) + cos(th(1) + th(2) + th(3)) * sin(th(0)) * sin(th(4))) + D(3) * cos(th(0)) - A(1) * cos(th(1)) * sin(th(0)) - D(4) * sin(th(1) + th(2) + th(3)) * sin(th(0)) - A(2) * cos(th(1)) * cos(th(2)) * sin(th(0)) + A(2) * sin(th(0)) * sin(th(1)) * sin(th(2)),
//         D(4) * (cos(th(4)) * sin(th(0)) - cos(th(1) + th(2) + th(3)) * cos(th(0)) * sin(th(4))) + D(3) * sin(th(0)) + A(1) * cos(th(0)) * cos(th(1)) + D(4) * sin(th(1) + th(2) + th(3)) * cos(th(0)) + A(2) * cos(th(0)) * cos(th(1)) * cos(th(2)) - A(2) * cos(th(0)) * sin(th(1)) * sin(th(2)),
//         0,
//         0,
//         0,
//         1;
//     Vector6d J2(6, 1); // second column of the jacobian matrix
//     J2 << -cos(th(0)) * (A(2) * sin(th(1) + th(2)) + A(1) * sin(th(1)) + D(4) * (sin(th(1) + th(2)) * sin(th(3)) - cos(th(1) + th(2)) * cos(th(3))) - D(4) * sin(th(4)) * (cos(th(1) + th(2)) * sin(th(3)) + sin(th(1) + th(2)) * cos(th(3)))),
//         -sin(th(0)) * (A(2) * sin(th(1) + th(2)) + A(1) * sin(th(1)) + D(4) * (sin(th(1) + th(2)) * sin(th(3)) - cos(th(1) + th(2)) * cos(th(3))) - D(4) * sin(th(4)) * (cos(th(1) + th(2)) * sin(th(3)) + sin(th(1) + th(2)) * cos(th(3)))),
//         A(2) * cos(th(1) + th(2)) - (D(4) * sin(th(1) + th(2) + th(3) + th(4))) / 2 + A(1) * cos(th(1)) + (D(4) * sin(th(1) + th(2) + th(3) - th(4))) / 2 + D(4) * sin(th(1) + th(2) + th(3)),
//         sin(th(0)),
//         -cos(th(0)),
//         0;
//     Vector6d J3(6, 1); // third column of the Jacobian matrix
//     J3 << cos(th(0)) * (D(4) * cos(th(1) + th(2) + th(3)) - A(2) * sin(th(1) + th(2)) + D(4) * sin(th(1) + th(2) + th(3)) * sin(th(4))),
//         sin(th(0)) * (D(4) * cos(th(1) + th(2) + th(3)) - A(2) * sin(th(1) + th(2)) + D(4) * sin(th(1) + th(2) + th(3)) * sin(th(4))),
//         A(2) * cos(th(1) + th(2)) - (D(4) * sin(th(1) + th(2) + th(3) + th(4))) / 2 + (D(4) * sin(th(1) + th(2) + th(3) - th(4))) / 2 + D(4) * sin(th(1) + th(2) + th(3)),
//         sin(th(0)),
//         -cos(th(0)),
//         0;
//     Vector6d J4(6, 1); // fourth column of the Jacobian matrix
//     J4 << D(4) * cos(th(0)) * (cos(th(1) + th(2) + th(3)) + sin(th(1) + th(2) + th(3)) * sin(th(4))),
//         D(4) * sin(th(0)) * (cos(th(1) + th(2) + th(3)) + sin(th(1) + th(2) + th(3)) * sin(th(4))),
//         D(4) * (sin(th(1) + th(2) + th(3) - th(4)) / 2 + sin(th(1) + th(2) + th(3)) - sin(th(1) + th(2) + th(3) + th(4)) / 2),
//         sin(th(0)),
//         -cos(th(0)),
//         0;
//     Vector6d J5(6, 1); // fifth column of the jacobian matrix
//     J5 << D(4) * cos(th(0)) * cos(th(1)) * cos(th(4)) * sin(th(2)) * sin(th(3)) - D(4) * cos(th(0)) * cos(th(1)) * cos(th(2)) * cos(th(3)) * cos(th(4)) - D(4) * sin(th(0)) * sin(th(4)) + D(4) * cos(th(0)) * cos(th(2)) * cos(th(4)) * sin(th(1)) * sin(th(3)) + D(4) * cos(th(0)) * cos(th(3)) * cos(th(4)) * sin(th(1)) * sin(th(2)),
//         D(4) * cos(th(0)) * sin(th(4)) + D(4) * cos(th(1)) * cos(th(4)) * sin(th(0)) * sin(th(2)) * sin(th(3)) + D(4) * cos(th(2)) * cos(th(4)) * sin(th(0)) * sin(th(1)) * sin(th(3)) + D(4) * cos(th(3)) * cos(th(4)) * sin(th(0)) * sin(th(1)) * sin(th(2)) - D(4) * cos(th(1)) * cos(th(2)) * cos(th(3)) * cos(th(4)) * sin(th(0)),
//         -D(4) * (sin(th(1) + th(2) + th(3) - th(4)) / 2 + sin(th(1) + th(2) + th(3) + th(4)) / 2),
//         sin(th(1) + th(2) + th(3)) * cos(th(0)),
//         sin(th(1) + th(2) + th(3)) * sin(th(0)),
//         -cos(th(1) + th(2) + th(3));
//     Vector6d J6(6, 1); // sixth column of the jacobian matrix
//     J6 << 0,
//         0,
//         0,
//         cos(th(4)) * sin(th(0)) - cos(th(1) + th(2) + th(3)) * cos(th(0)) * sin(th(4)),
//         -cos(th(0)) * cos(th(4)) - cos(th(1) + th(2) + th(3)) * sin(th(0)) * sin(th(4)),
//         -sin(th(1) + th(2) + th(3)) * sin(th(4));

//     // to obtain the jacobian matrix we just need to put in jacobianMatrix the above six vectors.
//     jacobianMatrix << J1, J2, J3, J4, J5, J6;
//     return jacobianMatrix;
// }

Eigen::Matrix4d directKin(Eigen::VectorXd Th)
{
    Matrix4d TransformationMatrix;
    // compute single transformations
    Matrix4d T10 = computeTransformationMatrix(Th(0), ALPHA(0), D(0), A(0));
    Matrix4d T21 = computeTransformationMatrix(Th(1), ALPHA(1), D(1), A(1));
    Matrix4d T32 = computeTransformationMatrix(Th(2), ALPHA(2), D(2), A(2));
    Matrix4d T43 = computeTransformationMatrix(Th(3), ALPHA(3), D(3), A(3));
    Matrix4d T54 = computeTransformationMatrix(Th(4), ALPHA(4), D(4), A(4));
    Matrix4d T65 = computeTransformationMatrix(Th(5), ALPHA(5), D(5), A(5));
    // compute transformation matrix from 0 to 6(T60=Tm)
    return TransformationMatrix = (T10) * (T21) * (T32) * (T43) * (T54) * (T65);
    // pe is the end effector postion: it's coincides with the first 3 rows of the last column( <3,1> vector) of Tm
    // pe = Tm.block<3, 1>(0, 3); // extract from row 0 and column 3 (so 4th column) a vector of 3x1 elements
    // Re is the rotation matrix-->3x3--> first 3 rows and first 3 columns of Tm
    // Re = Tm.block<3, 3>(0, 0);
}

Matrix4d computeTransformationMatrix(double th, double alpha, double distance, double a)
{
    Matrix4d Tij;
    Tij << cos(th), -sin(th) * cos(alpha), sin(th) * sin(alpha), a * cos(th),
        sin(th), cos(th) * cos(alpha), -cos(th) * sin(alpha), a * sin(th),
        0, sin(alpha), cos(alpha), distance,
        0, 0, 0, 1;
    return Tij;
}