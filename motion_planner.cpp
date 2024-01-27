#include "motion_planner.h"
using namespace Eigen;
using namespace std;

Matrix4d computeTransformationMatrix(double th, double alpha, double distance, double a)
{
    Matrix4d Tij;
    Tij << cos(th), -sin(th) * cos(alpha), sin(th) * sin(alpha), a * cos(th),
        sin(th), cos(th) * cos(alpha), -cos(th) * sin(alpha), a * sin(th),
        0, sin(alpha), cos(alpha), distance,
        0, 0, 0, 1;
    return Tij;
}

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

Matrix6d computeJacobian(Vector6d th)
{
    Matrix6d jacobianMatrix;
    Vector6d J1(6, 1); // first column of the jacobian matrix
    J1 << D(4) * (cos(th(0)) * cos(th(4)) + cos(th(1) + th(2) + th(3)) * sin(th(0)) * sin(th(4))) + D(3) * cos(th(0)) - A(1) * cos(th(1)) * sin(th(0)) - D(4) * sin(th(1) + th(2) + th(3)) * sin(th(0)) - A(2) * cos(th(1)) * cos(th(2)) * sin(th(0)) + A(2) * sin(th(0)) * sin(th(1)) * sin(th(2)),
        D(4) * (cos(th(4)) * sin(th(0)) - cos(th(1) + th(2) + th(3)) * cos(th(0)) * sin(th(4))) + D(3) * sin(th(0)) + A(1) * cos(th(0)) * cos(th(1)) + D(4) * sin(th(1) + th(2) + th(3)) * cos(th(0)) + A(2) * cos(th(0)) * cos(th(1)) * cos(th(2)) - A(2) * cos(th(0)) * sin(th(1)) * sin(th(2)),
        0,
        0,
        0,
        1;
    Vector6d J2(6, 1); // second column of the jacobian matrix
    J2 << -cos(th(0)) * (A(2) * sin(th(1) + th(2)) + A(1) * sin(th(1)) + D(4) * (sin(th(1) + th(2)) * sin(th(3)) - cos(th(1) + th(2)) * cos(th(3))) - D(4) * sin(th(4)) * (cos(th(1) + th(2)) * sin(th(3)) + sin(th(1) + th(2)) * cos(th(3)))),
        -sin(th(0)) * (A(2) * sin(th(1) + th(2)) + A(1) * sin(th(1)) + D(4) * (sin(th(1) + th(2)) * sin(th(3)) - cos(th(1) + th(2)) * cos(th(3))) - D(4) * sin(th(4)) * (cos(th(1) + th(2)) * sin(th(3)) + sin(th(1) + th(2)) * cos(th(3)))),
        A(2) * cos(th(1) + th(2)) - (D(4) * sin(th(1) + th(2) + th(3) + th(4))) / 2 + A(1) * cos(th(1)) + (D(4) * sin(th(1) + th(2) + th(3) - th(4))) / 2 + D(4) * sin(th(1) + th(2) + th(3)),
        sin(th(0)),
        -cos(th(0)),
        0;
    Vector6d J3(6, 1); // third column of the Jacobian matrix
    J3 << cos(th(0)) * (D(4) * cos(th(1) + th(2) + th(3)) - A(2) * sin(th(1) + th(2)) + D(4) * sin(th(1) + th(2) + th(3)) * sin(th(4))),
        sin(th(0)) * (D(4) * cos(th(1) + th(2) + th(3)) - A(2) * sin(th(1) + th(2)) + D(4) * sin(th(1) + th(2) + th(3)) * sin(th(4))),
        A(2) * cos(th(1) + th(2)) - (D(4) * sin(th(1) + th(2) + th(3) + th(4))) / 2 + (D(4) * sin(th(1) + th(2) + th(3) - th(4))) / 2 + D(4) * sin(th(1) + th(2) + th(3)),
        sin(th(0)),
        -cos(th(0)),
        0;
    Vector6d J4(6, 1); // fourth column of the Jacobian matrix
    J4 << D(4) * cos(th(0)) * (cos(th(1) + th(2) + th(3)) + sin(th(1) + th(2) + th(3)) * sin(th(4))),
        D(4) * sin(th(0)) * (cos(th(1) + th(2) + th(3)) + sin(th(1) + th(2) + th(3)) * sin(th(4))),
        D(4) * (sin(th(1) + th(2) + th(3) - th(4)) / 2 + sin(th(1) + th(2) + th(3)) - sin(th(1) + th(2) + th(3) + th(4)) / 2),
        sin(th(0)),
        -cos(th(0)),
        0;
    Vector6d J5(6, 1); // fifth column of the jacobian matrix
    J5 << D(4) * cos(th(0)) * cos(th(1)) * cos(th(4)) * sin(th(2)) * sin(th(3)) - D(4) * cos(th(0)) * cos(th(1)) * cos(th(2)) * cos(th(3)) * cos(th(4)) - D(4) * sin(th(0)) * sin(th(4)) + D(4) * cos(th(0)) * cos(th(2)) * cos(th(4)) * sin(th(1)) * sin(th(3)) + D(4) * cos(th(0)) * cos(th(3)) * cos(th(4)) * sin(th(1)) * sin(th(2)),
        D(4) * cos(th(0)) * sin(th(4)) + D(4) * cos(th(1)) * cos(th(4)) * sin(th(0)) * sin(th(2)) * sin(th(3)) + D(4) * cos(th(2)) * cos(th(4)) * sin(th(0)) * sin(th(1)) * sin(th(3)) + D(4) * cos(th(3)) * cos(th(4)) * sin(th(0)) * sin(th(1)) * sin(th(2)) - D(4) * cos(th(1)) * cos(th(2)) * cos(th(3)) * cos(th(4)) * sin(th(0)),
        -D(4) * (sin(th(1) + th(2) + th(3) - th(4)) / 2 + sin(th(1) + th(2) + th(3) + th(4)) / 2),
        sin(th(1) + th(2) + th(3)) * cos(th(0)),
        sin(th(1) + th(2) + th(3)) * sin(th(0)),
        -cos(th(1) + th(2) + th(3));
    Vector6d J6(6, 1); // sixth column of the jacobian matrix
    J6 << 0,
        0,
        0,
        cos(th(4)) * sin(th(0)) - cos(th(1) + th(2) + th(3)) * cos(th(0)) * sin(th(4)),
        -cos(th(0)) * cos(th(4)) - cos(th(1) + th(2) + th(3)) * sin(th(0)) * sin(th(4)),
        -sin(th(1) + th(2) + th(3)) * sin(th(4));

    // to obtain the jacobian matrix we just need to put in jacobianMatrix the above six vectors.
    jacobianMatrix << J1, J2, J3, J4, J5, J6;
    return jacobianMatrix;
}

Matrix4d directKin(Vector6d Th)
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

void inverseKin(MatrixXd &Th, double scaleFactor, VectorXd &pe60, Matrix3d &Re, Matrix4d &Tm)
{
    Tm.topLeftCorner<3, 3>() = Re;
    Tm.block<3, 1>(0, 3) = pe60;
    Tm.block<1, 4>(3, 0) << 0, 0, 0, 1;
    // find  theta1-->th1
    Vector4d p50 = Tm * Vector4d(0, 0, -D(5), 1); // p50 vector
    double psi = atan2(p50(1), p50(0));           // angle
    double p50norm = sqrt(p50(1) * p50(1) + p50(0) * p50(0));
    if (p50norm < D(3))
    {
        Th << numeric_limits<double>::quiet_NaN(),
            numeric_limits<double>::quiet_NaN(),
            numeric_limits<double>::quiet_NaN(),
            numeric_limits<double>::quiet_NaN(),
            numeric_limits<double>::quiet_NaN(),
            numeric_limits<double>::quiet_NaN();
        cout << "Position requested in the unreachable cylinder" << endl;
    }
    double phi1_1 = acos(D(3) / p50norm);
    double phi1_2 = -phi1_1;

    double th1_1 = psi + phi1_1 + M_PI_2;
    double th1_2 = psi + phi1_2 + M_PI_2;

    double p61z_1 = pe60(0) * sin(th1_1) - pe60(1) * cos(th1_1);
    double p61z_2 = pe60(0) * sin(th1_2) - pe60(1) * cos(th1_2);

    double th5_1_1 = acos((p61z_1 - D(3)) / D(5));
    double th5_1_2 = -acos((p61z_1 - D(3)) / D(5));
    double th5_2_1 = acos((p61z_2 - D(3)) / D(5));
    double th5_2_2 = -acos((p61z_2 - D(3)) / D(5));

    Matrix4d T10_1 = computeTransformationMatrix(th1_1, ALPHA(0), D(0), A(0));
    Matrix4d T10_2 = computeTransformationMatrix(th1_2, ALPHA(0), D(0), A(0));

    Matrix4d T16_1 = ((T10_1.inverse()) * Tm).inverse();
    Matrix4d T16_2 = ((T10_2.inverse()) * Tm).inverse();

    double zy_1 = T16_1(1, 2);
    double zx_1 = T16_1(0, 2);

    double zy_2 = T16_2(1, 2);
    double zx_2 = T16_2(0, 2);
    double th6_1_1;
    double th6_1_2;
    double th6_2_1;
    double th6_2_2;

    if (almostZero(sin(th5_1_1)) || (almostZero(zy_1) && almostZero(zx_1)))
    {
        cout << "singular configuration. Choosing arbitrary th6" << endl;
        th6_1_1 = 0;
    }
    else
    {
        th6_1_1 = atan2((-zy_1 / sin(th5_1_1)), (zx_1 / sin(th5_1_1)));
    }

    if (almostZero(sin(th5_1_2)) || (almostZero(zy_1) && almostZero(zx_1)))
    {
        cout << "singular configuration. Choosing arbitrary th6" << endl;
        th6_1_2 = 0;
    }
    else
    {
        th6_1_2 = atan2((-zy_1 / sin(th5_1_2)), (zx_1 / sin(th5_1_2)));
    }

    if (almostZero(sin(th5_2_1)) || (almostZero(zy_2) && almostZero(zx_2)))
    {
        cout << "singular configuration. Choosing arbitrary th6" << endl;
        th6_2_1 = 0;
    }
    else
    {
        th6_2_1 = atan2((-zy_2 / sin(th5_2_1)), (zx_2 / sin(th5_2_1)));
    }

    if (almostZero(sin(th5_2_2)) || (almostZero(zy_2) && almostZero(zx_2)))
    {
        cout << "singular configuration. Choosing arbitrary th6" << endl;
        th6_2_2 = 0;
    }
    else
    {
        th6_2_2 = atan2((-zy_2 / sin(th5_2_2)), (zx_2 / sin(th5_2_2)));
    }

    Matrix4d T61_1 = T16_1.inverse();
    Matrix4d T61_2 = T16_2.inverse();

    Matrix4d T54_1_1 = computeTransformationMatrix(th5_1_1, ALPHA(4), D(4), A(4));
    Matrix4d T54_1_2 = computeTransformationMatrix(th5_1_2, ALPHA(4), D(4), A(4));
    Matrix4d T54_2_1 = computeTransformationMatrix(th5_2_1, ALPHA(4), D(4), A(4));
    Matrix4d T54_2_2 = computeTransformationMatrix(th5_2_2, ALPHA(4), D(4), A(4));

    Matrix4d T65_1_1 = computeTransformationMatrix(th6_1_1, ALPHA(5), D(5), A(5));
    Matrix4d T65_1_2 = computeTransformationMatrix(th6_1_2, ALPHA(5), D(5), A(5));
    Matrix4d T65_2_1 = computeTransformationMatrix(th6_2_1, ALPHA(5), D(5), A(5));
    Matrix4d T65_2_2 = computeTransformationMatrix(th6_2_2, ALPHA(5), D(5), A(5));

    Matrix4d T41_1_1 = T61_1 * (T54_1_1 * T65_1_1).inverse();
    Matrix4d T41_1_2 = T61_1 * (T54_1_2 * T65_1_2).inverse();
    Matrix4d T41_2_1 = T61_2 * (T54_2_1 * T65_2_1).inverse();
    Matrix4d T41_2_2 = T61_2 * (T54_2_2 * T65_2_2).inverse();

    // Point position depending on 4 angles found previously
    Vector3d P31_1_1, P31_1_2, P31_2_1, P31_2_2;
    Vector4d P;
    P = T41_1_1 * Vector4d(0, -D(3), 0, 1);
    P31_1_1 = P.head(3);
    P = T41_1_2 * Vector4d(0, -D(3), 0, 1);
    P31_1_2 = P.head(3);
    P = T41_2_1 * Vector4d(0, -D(3), 0, 1);
    P31_2_1 = P.head(3);
    P = T41_2_2 * Vector4d(0, -D(3), 0, 1);
    P31_2_2 = P.head(3);

    double C;

    double th3_1_1_1;
    double th3_1_1_2;
    C = (((P31_1_1).norm() * (P31_1_1).norm()) - A(2) * A(2) - A(3) * A(3)) / (2 * A(2) * A(3));
    if (abs(C) > 1)
    {
        cout << "Point out of the work space" << endl;
        th3_1_1_1 = numeric_limits<double>::quiet_NaN();
        th3_1_1_2 = numeric_limits<double>::quiet_NaN();
    }
    else
    {
        th3_1_1_1 = acos(C);
        th3_1_1_2 = -acos(C);
    }

    double th3_1_2_1;
    double th3_1_2_2;
    C = (((P31_1_2).norm() * (P31_1_1).norm()) - A(2) * A(2) - A(3) * A(3)) / (2 * A(2) * A(3));
    if (abs(C) > 1)
    {
        cout << "Point out of the work space" << endl;
        th3_1_2_1 = numeric_limits<double>::quiet_NaN();
        th3_1_2_2 = numeric_limits<double>::quiet_NaN();
    }
    else
    {
        th3_1_2_1 = acos(C);
        th3_1_2_2 = -acos(C);
    }

    double th3_2_1_1;
    double th3_2_1_2;
    C = (((P31_2_1).norm() * (P31_1_1).norm()) - A(2) * A(2) - A(3) * A(3)) / (2 * A(2) * A(3));
    if (abs(C) > 1)
    {
        cout << "Point out of the work space" << endl;
        th3_2_1_1 = numeric_limits<double>::quiet_NaN();
        th3_2_1_2 = numeric_limits<double>::quiet_NaN();
    }
    else
    {
        th3_2_1_1 = acos(C);
        th3_2_1_2 = -acos(C);
    }

    double th3_2_2_1;
    double th3_2_2_2;
    C = (((P31_2_2).norm() * (P31_1_1).norm()) - A(2) * A(2) - A(3) * A(3)) / (2 * A(2) * A(3));
    if (abs(C) > 1)
    {
        cout << "Point out of the work space" << endl;
        th3_2_2_1 = numeric_limits<double>::quiet_NaN();
        th3_2_2_2 = numeric_limits<double>::quiet_NaN();
    }
    else
    {
        th3_2_2_1 = acos(C);
        th3_2_2_2 = -acos(C);
    }

    double th2_1_1_1 = -atan2(P31_1_1(1), -P31_1_1(0)) + asin((A(2) * sin(th3_1_1_1)) / P31_1_1.norm());
    double th2_1_1_2 = -atan2(P31_1_1(1), -P31_1_1(0)) + asin((A(2) * sin(th3_1_1_2)) / P31_1_1.norm());
    double th2_1_2_1 = -atan2(P31_1_2(1), -P31_1_2(0)) + asin((A(2) * sin(th3_1_2_1)) / P31_1_2.norm());
    double th2_1_2_2 = -atan2(P31_1_2(1), -P31_1_2(0)) + asin((A(2) * sin(th3_1_2_2)) / P31_1_2.norm());
    double th2_2_1_1 = -atan2(P31_2_1(1), -P31_2_1(0)) + asin((A(2) * sin(th3_2_1_1)) / P31_2_1.norm());
    double th2_2_1_2 = -atan2(P31_2_1(1), -P31_2_1(0)) + asin((A(2) * sin(th3_2_1_2)) / P31_2_1.norm());
    double th2_2_2_1 = -atan2(P31_2_2(1), -P31_2_2(0)) + asin((A(2) * sin(th3_2_2_1)) / P31_2_2.norm());
    double th2_2_2_2 = -atan2(P31_2_2(1), -P31_2_2(0)) + asin((A(2) * sin(th3_2_2_2)) / P31_2_2.norm());

    Matrix4d T21, T32, T41, T43;
    double xy, xx;

    T21 = computeTransformationMatrix(th2_1_1_1, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_1_1_1, ALPHA(2), D(2), A(2));
    T41 = T41_1_1;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_1_1_1 = atan2(xy, xx);

    T21 = computeTransformationMatrix(th2_1_1_2, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_1_1_2, ALPHA(2), D(2), A(2));
    T41 = T41_1_1;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_1_1_2 = atan2(xy, xx);

    T21 = computeTransformationMatrix(th2_1_2_1, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_1_2_1, ALPHA(2), D(2), A(2));
    T41 = T41_1_2;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_1_2_1 = atan2(xy, xx);

    T21 = computeTransformationMatrix(th2_1_2_2, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_1_2_2, ALPHA(2), D(2), A(2));
    T41 = T41_1_2;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_1_2_2 = atan2(xy, xx);

    T21 = computeTransformationMatrix(th2_2_1_1, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_2_1_1, ALPHA(2), D(2), A(2));
    T41 = T41_2_1;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_2_1_1 = atan2(xy, xx);

    T21 = computeTransformationMatrix(th2_2_1_2, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_2_1_2, ALPHA(2), D(2), A(2));
    T41 = T41_2_1;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_2_1_2 = atan2(xy, xx);

    T21 = computeTransformationMatrix(th2_2_2_1, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_2_2_1, ALPHA(2), D(2), A(2));
    T41 = T41_2_2;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_2_2_1 = atan2(xy, xx);

    T21 = computeTransformationMatrix(th2_2_2_2, ALPHA(1), D(1), A(1));
    T32 = computeTransformationMatrix(th3_2_2_2, ALPHA(2), D(2), A(2));
    T41 = T41_2_2;
    T43 = ((T21 * T32).inverse()) * T41;
    xy = T43(1, 0);
    xx = T43(0, 0);
    double th4_2_2_2 = atan2(xy, xx);
    cout << "hi1" << endl;
    // Finally compute the inverse kinematic and find the joint angles (remembering we started from final rotation matrix)
    Th << th1_1, th1_1, th1_1, th1_1, th1_2, th1_2, th1_2, th1_2,
        th2_1_1_1, th2_1_1_2, th2_1_2_1, th2_1_2_2, th2_2_2_1, th2_2_1_2, th2_2_2_1, th2_2_2_2,
        th3_1_1_1, th3_1_1_2, th3_1_2_1, th3_1_2_2, th3_2_2_1, th3_2_1_2, th3_2_2_1, th3_2_2_2,
        th4_1_1_1, th4_1_1_2, th4_1_2_1, th4_1_2_2, th4_2_2_1, th4_2_1_2, th4_2_2_1, th4_2_2_2,
        th5_1_1, th5_1_1, th5_1_2, th5_1_2, th5_2_2, th5_2_1, th5_2_2, th5_2_2,
        th6_1_1, th6_1_1, th6_1_2, th6_1_2, th6_2_2, th6_2_1, th6_2_2, th6_2_2;
}

Vector3d p_i(double time, Vector3d p1, Vector3d p2)
{
    // since interpolation evaluates which weight must be given to each vectors, we need to normalize the time with reference with the time delta of the path (path_dt)
    const double normalizedTime = time / path_dt;
    if (normalizedTime < 1)
        return (normalizedTime * p2) + ((1 - normalizedTime) * p1);
    else
        return p2; // zero weight is given to p1, since time exceed our delta
}

Quaterniond q_i(double time, Quaterniond q1, Quaterniond q2)
{
    const double normalizedTime = time / path_dt;
    if (normalizedTime < 1)
        return (q1.slerp(normalizedTime, q2));
    else
        return q2;
}

Path inverseDiffQuat(Vector8d UR5measures, Vector3d p_start, Vector3d p_end, Quaterniond q_start, Quaterniond q_end)
{
    // gripper state
    Vector2d gripper_state{UR5measures(6), UR5measures(7)}; // last two rows of measures vector
    // joint values at time=t, their derivatives in time and (velocities+error(e))-->real velocities
    Vector6d jointValues, jointValues_dot, correctedVelocities;
    Path pathToFollow; //-->returned path
    Matrix4d transfromationMatrixAtT;
    Vector3d positionAtT, linVelocitiesAtT, angVelocitiesAtT;
    Quaterniond quaternionAtT, quaternionAngVelAtT, quaternionErr;
    Matrix6d jacobianAtT, pseudoRightInverseJacAtT; // jacobian and it pseudo right inverse matrices

    Matrix3d Kp, Kq;         // are the matrices used to correction of position and quaternion
    Kp = Kp.Identity() * 10; // 10 is the correction factor
    Kq = Kq.Identity() * 1;  // 1 is the correction factor

    // add point to path
    for (int i = 0; i < 6; ++i)
    {
        jointValues(i) = UR5measures(i);
    }

    // at each dt compute new joints state
    for (double t = dt; t < path_dt; t += dt)
    {
        // compute direct kinematics a time=T
        transfromationMatrixAtT = directKin(jointValues);
        positionAtT = transfromationMatrixAtT.block(0, 3, 3, 1);
        quaternionAtT = Quaterniond(transfromationMatrixAtT.block(0, 0, 3, 3));

        // compute velocities at time=T
        linVelocitiesAtT = (p_i(t, p_start, p_end) - p_i(t - dt, p_start, p_end)) / dt;
        quaternionAngVelAtT = q_i(t + dt, q_start, q_end) * q_i(t, q_start, q_end).conjugate();
        angVelocitiesAtT = (quaternionAngVelAtT.vec() * 2) / dt;

        jacobianAtT = computeJacobian(jacobianAtT);
        pseudoRightInverseJacAtT = ((jacobianAtT.transpose() * jacobianAtT + jacobianAtT.Identity() * 0.0001).inverse()) * jacobianAtT.transpose();

        if (abs(jacobianAtT.determinant()) < 0.00001)
        {
            ROS_WARN("Near singular configuration");
        }

        quaternionErr = q_i(t, q_start, q_end) * quaternionAtT.conjugate();
        Vector3d linearPosError = p_i(t, p_start, p_end) - positionAtT;

        correctedVelocities << linVelocitiesAtT + (Kp * linearPosError), angVelocitiesAtT + (Kq * quaternionErr.vec());

        jointValues_dot = pseudoRightInverseJacAtT * correctedVelocities;
        jointValues = jointValues + (jointValues_dot * dt);

        pathToFollow = insert_new_path_instance(pathToFollow, jointValues, gripper_state);
    }
    return pathToFollow;
}

Path insert_new_path_instance(Path p, Vector6d js, Vector2d gripper_state)
    {
        p.conservativeResize(p.rows() + 1, p.cols());
        p.row(p.rows() - 1) = V8d{js(0), js(1), js(2), js(3), js(4), js(5), gripper_state(0), gripper_state(1)};
        return p;
    }

