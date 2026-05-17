#include "caams.hpp"
#include <math.h>

namespace caams
{
	Eigen::Matrix3x4d G(Eigen::Vector4d const & p){
		Eigen::Matrix3x4d r;
		r << -p(1), p(0),-p(3), p(2),
			 -p(2), p(3), p(0),-p(1),
			 -p(3),-p(2), p(1), p(0);
		return r;
	}
	
	Eigen::Matrix3x4d L(Eigen::Vector4d const & p){
		Eigen::Matrix3x4d r;
		r << -p(1), p(0), p(3),-p(2),
			 -p(2),-p(3), p(0), p(1),
			 -p(3), p(2),-p(1), p(0);
		return r;
	}
	
	Eigen::Matrix3d SS(Eigen::Vector3d const & a){
		Eigen::Matrix3d r;
		r << 0.0 ,-a(2), a(1),
			 a(2), 0.0 ,-a(0),
			-a(1), a(0), 0.0;
		return r;
	}
	
	Eigen::Vector4d pAA(double angle, Eigen::Vector3d const & axis){
		Eigen::Vector4d r;
		angle/=2.0;
		r(0) = cos(angle);
		r.segment<3>(1) = axis.normalized()*sin(angle);
		return r;
	}

    Eigen::Vector4d pA(Eigen::Matrix3d const & A){
        double trace = A(0,0) + A(1,1) + A(2,2);
        double e0 = (trace+1)/4.0;
        double e1;
        double e2;
        double e3;
        if(e0<0.0)e0=0.0;
        e0 = sqrt(e0);
        if(e0<1e-6){
            e0 = 0.0;
            double e1s = (1 + 2.0*A(0,0) - trace)/4.0;
            double e2s = (1 + 2.0*A(1,1) - trace)/4.0;
            double e3s = (1 + 2.0*A(2,2) - trace)/4.0;
            int i_max;
            if(e1s>e2s){
                if(e1s>e3s){
                    i_max = 1;
                }else{
                    i_max = 3;
                }
            }else{
                if(e2s>e3s){
                    i_max = 2;
                }else{
                    i_max = 3;
                }
            }
            switch(i_max){
            case 1:
                e1 = sqrt(e1s);
                e2 = (A(1,0)+A(0,1))/(4.0*e1);
                e3 = (A(2,0)+A(0,2))/(4.0*e1);
                break;
            case 2:
                e2 = sqrt(e2s);
                e1 = (A(1,0)+A(0,1))/(4.0*e2);
                e3 = (A(2,1)+A(1,2))/(4.0*e2);
                break;
            case 3:
                e3 = sqrt(e3s);
                e1 = (A(2,0)+A(0,2))/(4.0*e3);
                e2 = (A(2,1)+A(1,2))/(4.0*e3);
                break;
            }
        }else{
            e1 = (A(2,1) - A(1,2))/(4.0*e0);
            e2 = (A(0,2) - A(2,0))/(4.0*e0);
            e3 = (A(1,0) - A(0,1))/(4.0*e0);
        }
        Eigen::Vector4d r;
        r << e0, e1, e2, e3;
        return r;
    }

    Eigen::Matrix3d Ap(Eigen::Vector4d const & p){
        //return G(p)*L(p).transpose();
        double e0_2 = p(0)*p(0);
        double e1_2 = p(1)*p(1);
        double e2_2 = p(2)*p(2);
        double e3_2 = p(3)*p(3);
        double t_e0_e1 = 2.0*p(0)*p(1);
        double t_e0_e2 = 2.0*p(0)*p(2);
        double t_e0_e3 = 2.0*p(0)*p(3);
        double t_e1_e2 = 2.0*p(1)*p(2);
        double t_e1_e3 = 2.0*p(1)*p(3);
        double t_e2_e3 = 2.0*p(2)*p(3);
        Eigen::Matrix3d A;
        A(0,0) = e0_2 + e1_2 - e2_2 - e3_2;
        A(1,0) = t_e0_e3 + t_e1_e2;
        A(2,0) = t_e1_e3 - t_e0_e2;
        A(0,1) = t_e1_e2 - t_e0_e3;
        A(1,1) = e0_2 - e1_2 + e2_2 - e3_2;
        A(2,1) = t_e2_e3 + t_e0_e1;
        A(0,2) = t_e1_e3 + t_e0_e2;
        A(1,2) = t_e2_e3 - t_e0_e1;
        A(2,2) = e0_2 - e1_2 - e2_2 + e3_2;
        return A;
    }

    Eigen::Matrix3d AAA(double angle, Eigen::Vector3d const & axis)
    {
        return Ap(pAA(angle, axis));
    }

    Eigen::Vector4d p_dot_omega(Eigen::Vector4d const & p, Eigen::Vector3d const & omega)
    {
        Eigen::Vector4d p_dot = 0.5*G(p).transpose()*omega;
        return p_dot;
    }

    Eigen::Vector4d p_dot_omega_p(Eigen::Vector4d const & p, Eigen::Vector3d const & omega_p)
    {
        Eigen::Vector4d p_dot = 0.5*L(p).transpose()*omega_p;
        return p_dot;
    }

    Eigen::Vector3d omega_p_dot(Eigen::Vector4d const & p, Eigen::Vector4d const & p_dot)
    {
        Eigen::Vector3d omega = 2.0*G(p)*p_dot;
        return omega;
    }

    Eigen::Vector3d omega_p_p_dot(Eigen::Vector4d const & p, Eigen::Vector4d const & p_dot)
    {
        Eigen::Vector3d omega_p = 2.0*L(p)*p_dot;
        return omega_p;
    }

    Eigen::Vector4d p_ddot_omega_dot(Eigen::Vector4d const & p, Eigen::Vector3d const & omega_dot, Eigen::Vector3d const & omega)
    {
        double w2 = omega.squaredNorm();
        return 0.5*G(p).transpose()*omega_dot - 0.25*w2*p;
    }

    Eigen::Vector4d p_ddot_omega_p_dot(Eigen::Vector4d const & p, Eigen::Vector3d const & omega_p_dot, Eigen::Vector3d const & omega_p)
    {
        double w2 = omega_p.squaredNorm();
        return 0.5*L(p).transpose()*omega_p_dot - 0.25*w2*p;
    }

    Eigen::Vector4d p_ddot_solve(
            Eigen::Vector4d &p,
            Eigen::Vector4d &p_dot,
            Eigen::Matrix3d &Jp,
            Eigen::Vector3d &np)
    {
        Eigen::Matrix4d H = 4.0*caams::L(p_dot).transpose()*Jp*caams::L(p);
        Eigen::Matrix4d R;
        R.topLeftCorner(3,4) = caams::L(p)*H;
        R.row(3) = p_dot.transpose();

        Eigen::Vector4d rhs;
        rhs << np, 0.0;

        Eigen::Vector4d y = rhs - R*p_dot;

        Eigen::Matrix4d A;
        A.topLeftCorner(3,4) = 2.0*Jp*caams::L(p);
        A.row(3) = p.transpose();

        Eigen::Vector4d x = A.inverse()*y;

        return x;
    }

    Eigen::Matrix4d a_plus(Eigen::Vector3d const & a){
		Eigen::Matrix4d r;
		r <<
			 0.0 ,-a(0),-a(1),-a(2),
			 a(0), 0.0 ,-a(2), a(1),
			 a(1), a(2), 0.0 ,-a(0),
			 a(2),-a(1), a(0), 0.0;
		return r;
	}

	Eigen::Matrix4d a_minus(Eigen::Vector3d const & a){
		Eigen::Matrix4d r;
		r <<
			 0.0 ,-a(0),-a(1),-a(2),
			 a(0), 0.0 , a(2),-a(1),
			 a(1),-a(2), 0.0 , a(0),
			 a(2), a(1),-a(0), 0.0;
		return r;
	}

	Eigen::Matrix3d J_p_cylinder_x_axis(double m, double r, double l){
        r*=r;
        l*=l;
        double Jxx = m*r/2.0;
        double Jyy = m/12.0*(3*r+l);
        double Jzz = Jyy;
		Eigen::Matrix3d res;
		res<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return res;
    }

	Eigen::Matrix3d J_p_cylinder_y_axis(double m, double r, double h){
        r*=r;
        h*=h;
        double Jxx = m/12.0*(3*r+h);
        double Jyy = m*r/2.0;
        double Jzz = Jxx;
		Eigen::Matrix3d res;
		res<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return res;
	}

	Eigen::Matrix3d J_p_cylinder_z_axis(double m, double r, double h){
        r*=r;
        h*=h;
        double Jxx = m/12.0*(3*r+h);
        double Jzz = m*r/2.0;
        double Jyy = Jxx;
		Eigen::Matrix3d res;
		res<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return res;
	}

	Eigen::Matrix3d J_p_sphere(double mass, double radius){
        double I_sphere = 2.0/5.0*mass*radius*radius;
		Eigen::Matrix3d r;
		r<<
			I_sphere, 0.0, 0.0,
			0.0, I_sphere, 0.0,
			0.0, 0.0, I_sphere;
		return r;
	}

	Eigen::Matrix3d J_p_cuboid(double mass, double dx, double dy, double dz){
        dx*=dx;
        dy*=dy;
        dz*=dz;
        mass*=1.0/12.0;
        double Jxx = mass*(dy+dz);
        double Jyy = mass*(dx+dz);
        double Jzz = mass*(dx+dy);
		Eigen::Matrix3d r;
		r<<
			Jxx, 0.0, 0.0,
			0.0, Jyy, 0.0,
			0.0, 0.0, Jzz;
		return r;
	}
}
