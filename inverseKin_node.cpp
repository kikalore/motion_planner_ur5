#include "inverseKin_node.h"
using namespace Eigen;
using namespace std;
/*INVERSE KINEMATICS COMPUTATION OF UR5
We have defined the followed:
->Th: the ur5 six joint angles to FIND
->pe: desidered cartesian position of the end effector
->Re: desidered rotation matrix of the end effector
*/
int main()
{
    // position and orientation desidered of the end-effector
    VectorXd pe60(3);
    pe60 << 1, 1, 1; // to modify with the desidered position-->vision node
    Matrix3d Re;
    Re << 1, 1, 1,
        1, 1, 1,
        1, 1, 1; // to modify with the desidered orientation-->vision node

    // joint angles to find
    MatrixXd JointAngles(6,8);

    // final transformation matrix
    Matrix4d Tm;

    // compute inverse kinematics
    inverseKin(JointAngles,1.0, pe60, Re, Tm);

    cout << "Possible configurations of joint angles:\n"<< JointAngles << endl; 
    return 0;
}