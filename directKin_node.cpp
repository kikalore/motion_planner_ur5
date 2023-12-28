#include "directKin_node.h"
using namespace Eigen;
using namespace std;
/*DIRECT KINEMATICS COMPUTATION OF UR5
We have defined the followed:
->Th: the ur5 six joint angles
->pe: cartesian position of the end effector
->Re: rotation matrix of the end effector
*/

int main()
{
    VectorXd Th(6), pe(3);
    Matrix3d Re;
    Matrix4d Tm;
    Th<<0,0,0,0,0,0;
    double scaleFactor=10;
    //direct kinematics
    directKin(Th,scaleFactor,pe,Re,Tm);
    // Print direct kinematics matrix
    print_eigen("DirectKinematics", Tm);
    cout<<endl; 
    return 0;
}