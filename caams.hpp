#ifndef CAAMS_H
#define CAAMS_H

#include <Eigen/Dense>

namespace Eigen
{
	typedef Matrix<double, 3, 4> Matrix3x4d;
}

namespace caams
{
	Eigen::Matrix3x4d G(Eigen::Vector4d const & p);
	Eigen::Matrix3x4d L(Eigen::Vector4d const & p);
	Eigen::Matrix3d SS(Eigen::Vector3d const & a);
	Eigen::Vector4d pAA(double angle, Eigen::Vector3d const & axis);
	Eigen::Matrix3d Ap(Eigen::Vector4d const & p);
    Eigen::Vector4d pA(Eigen::Matrix3d const & A);
    Eigen::Matrix3d AAA(double angle, Eigen::Vector3d const & axis);
    Eigen::Matrix4d a_plus(Eigen::Vector3d const & a);
	Eigen::Matrix4d a_minus(Eigen::Vector3d const & a);
    Eigen::Vector4d p_dot_omega(Eigen::Vector4d const & p, Eigen::Vector3d const & omega);
    Eigen::Vector4d p_dot_omega_p(Eigen::Vector4d const & p, Eigen::Vector3d const & omega_p);
    Eigen::Vector3d omega_p_dot(Eigen::Vector4d const & p, Eigen::Vector4d const & p_dot);
    Eigen::Vector3d omega_p_p_dot(Eigen::Vector4d const & p, Eigen::Vector4d const & p_dot);
    Eigen::Vector4d p_ddot_omega_dot(Eigen::Vector4d const & p, Eigen::Vector3d const & omega_dot, Eigen::Vector3d const & omega);
    Eigen::Vector4d p_ddot_omega_p_dot(Eigen::Vector4d const & p, Eigen::Vector3d const & omega_p_dot, Eigen::Vector3d const & omega_p);
    Eigen::Vector4d p_ddot_solve(
            Eigen::Vector4d &p,
            Eigen::Vector4d &p_dot,
            Eigen::Matrix3d &Jp,
            Eigen::Vector3d &np);

    // some inertial tensor functions
	Eigen::Matrix3d J_p_cylinder_x_axis(double m, double r, double l);
	Eigen::Matrix3d J_p_cylinder_y_axis(double m, double r, double h);
	Eigen::Matrix3d J_p_cylinder_z_axis(double m, double r, double h);
	Eigen::Matrix3d J_p_sphere(double mass, double radius);
	Eigen::Matrix3d J_p_cuboid(double mass, double dx, double dy, double dz);

}

#endif
