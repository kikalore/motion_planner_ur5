#include "motion_planner.h"
using namespace Eigen;
using namespace std;
// testing
int main()
{
    Matrix6d myJac;
    Matrix4d directMatrix, directMatrix2;
    Vector6d myJointVariables;
    Vector6d myvars;
    Vector3d pe60;
    Matrix3d Re;
    Matrix4d Tm;
    // testing with homing position variables
    myJointVariables << -0.0001, -0.00016, -0.00006, -0.00014, -0.00004, 0.00001;
    myvars<< 1,1,1,1,1,1;

    // matrix from base to world
    cout << "\nrototransMatrix\n"
         << rotation_180z_axis_and_offset() << endl;

    // testing direct kinematics
    directMatrix = directKin(myvars);
    print_eigen("Direct matrix with reference to BASE frame\n", directMatrix);
    cout << endl;

    // PREMULTIPLY to express position in a FIXED frame
    directMatrix = rotation_180z_axis_and_offset() * directKin(myvars);
    print_eigen("Direct matrix with reference to WORLD frame (NOT CORRECT?)\n", directMatrix);
    cout << endl;

    // testing Jacobian computation
    myJac = computeJacobian(myvars);
    print_eigen("Jacobian matrix\n", myJac);

    // testing interpolation (position and quaternion)
    Quaterniond q1(8, 2, 3, 4);
    Quaterniond q2(5, 6, 7, 8);
    Vector3d p1(1, 2, 3);
    Vector3d p2(4, 5, 6);
    // test case 1: tempo normalizzato minore di 1
    double time1 = 5; // Tempo arbitrario-->normalized time will be 0.5, since it's < 1 we expected slerp as a result
    Quaterniond result1 = myslerp(time1, q1, q2);
    Vector3d result_p1 = lerp(time1, p1, p2);
    cout << "Test case 1:" << endl;
    cout << "Result: " << result1.coeffs().transpose() << "\t" << result_p1.transpose() << endl;
    cout << endl;
    // test case 2: tempo normalizzato maggiore di 1
    double time2 = 30; // Tempo arbitrario-->normalized time will be 3, since it's > 1 we expected q2 as a result
    Quaterniond result2 = myslerp(time2, q1, q2);
    Quaterniond expected2 = q2;
    Vector3d result_p2 = lerp(time2, p1, p2);
    Vector3d expected_p2 = p2;
    cout << "Test case 2:" << endl;
    cout << "Result: " << result2.coeffs().transpose() << "\t" << result_p2.transpose() << endl;
    cout << endl;

    // testing inverse kinematics
    Matrix<double, 6, 8> Th;
    Vector3d my_pos60 = directMatrix.col(3).head(3);
    Matrix3d my_R60 = directMatrix.block<3, 3>(0, 0);
    Th = inverseKin(my_pos60, my_R60, 10);
    print_eigen("joint variables from inverse kinematics\n", Th);
    return 0;
}
