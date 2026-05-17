#include "conversions.h"
#include "caams.hpp"
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <iostream>
#include <vector>

using namespace boost::numeric::odeint;
using namespace Eigen;
typedef std::vector<double> state_type;
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;

class Ballistics
{
public:
    double m;
    double v0;
    double vw;
    double heading_range;
    double heading_wind;
    double diameter;
    double length;
    double d_cp;
    double d_rr;
    double radius;
    double A;
    double g;
    Vector3d a_g_vec;
    double h_sight;
    Matrix3d Jp;
    Vector3d r_cp;
    double k_rr;

    double R;
    double Md;
    double Mv;
    double L;
    double p0;
    double Tc;
    double T0;
    double humidity;
    double rho;
    Vector3d vw_vec;
    MatrixXd Cd_table;
    double c;

    double psat(double T)
    {
        T -= 273.15;
        return 611.21*exp((18.678-T/234.5)*(T/(257.14+T)));
    }

    double rho_from_weather(double T, double p, double humid)
    {
        double pv = humid/100.0*psat(T);
        double pd = p - pv;
        return (pd*Md + pv*Mv)/R/T;
    }

    Vector3d wind_vec(double heading_range, double heading_wind, double v_wind)
    {
        double theta = heading_wind - heading_range;
        theta *= M_PI/180.0;
        Vector3d r{cos(theta)*v_wind, 0.0, sin(theta)*v_wind};
        return r;
    }

    double Cd_at_mach(double m)
    {
        if(m<0.0 || m>5.0){
            std::cout << "Cd_at_mach: m out of range m=" << m << std::endl;
            throw;
        }
        //std::cout << "m:" << m << std::endl;
        for(int i=0;;i++){
            if( m>=Cd_table(0,i) && m<=Cd_table(0,i+1)){
                double alpha = (m-Cd_table(0,i))/(Cd_table(0,i+1)-Cd_table(0,i));
                return Cd_table(1,i)*(1.0-alpha) + Cd_table(1,i+1)*alpha;
            }
        }
    }

    double speed_of_sound(double T)
    {
        return 20.047*sqrt(T);
    }

    Ballistics(void)
    {
        m = 95*kg_p_gr;
        v0 = 3300*mps_p_fps;
        vw = 0*mps_p_kmph;
        heading_range = 0.0;
        heading_wind = 135.0;
        diameter = 6.7e-3;
        length = 50.0e-3;
        d_cp = 10.0e-3;
        d_rr = 8*m_p_inch;
        radius = diameter/2;
        A = M_PI*radius*radius;
        g = 9.81;
        a_g_vec << 0.0, -g, 0.0;
        h_sight = 1.5*m_p_inch;
        Jp = caams::J_p_cylinder_x_axis(m, radius, length);
        r_cp << d_cp, 0.0, 0.0;
        k_rr = 2*M_PI/d_rr;

        R = 8.31466;
        Md = 0.0289652;
        Mv = 0.018016;
        L = 0.0065;
        p0 = 29.53*Pa_p_inHg;
        Tc = 15.0;
        T0 = Tc + 273.15;
        humidity = 78;
        rho = rho_from_weather(T0, p0, humidity);
        vw_vec = wind_vec(heading_range, heading_wind, vw);
        Cd_table.resize(2,9);
        Cd_table << 0.,   0.7,  0.8,  0.9,  1.0, 1.1,   1.6,  3.5,  5.0,
                    0.12, 0.12, 0.13, 0.15, 0.4, 0.405, 0.33, 0.22, 0.17;
        c = speed_of_sound(T0);

    }

    //
    // This class is a functor and the function computes the derivatives
    // of the state vector.
    //
    void operator() (const state_type &y_in, state_type &dy_in, const double t)
    {
        Map<const VectorXd> y(y_in.data(), 14);
        Map<VectorXd> dy(dy_in.data(), 14);
        Vector3d v_vec = y.segment<3>(3);
        Vector4d p = y.segment<4>(6);
        Vector4d pdot = y.segment<4>(10);
        p.normalize();
        Matrix3d A_g = caams::Ap(p);
        Vector3d v_rel = v_vec + vw_vec;
        double v = v_rel.norm();
        double mach = v/c;
        double c_d = Cd_at_mach(mach);
        Vector3d Cd_factor{1.0, 2.0, 2.0};
        Vector3d v_air_obj = -(A_g.transpose()*v_rel);
        Vector3d n_air_obj = v_air_obj/v;
        Vector3d F_air_obj = 0.5*rho*v*v*c_d*A*n_air_obj;
        F_air_obj.array() *= Cd_factor.array();
        Vector3d F_air = A_g * F_air_obj;
        Vector3d a_vec = a_g_vec + F_air/m;
        Vector3d np = r_cp.cross(F_air_obj);
        Vector4d pddot = caams::p_ddot_solve(p,pdot,Jp,np);

        dy.segment<3>(0) = v_vec;
        dy.segment<3>(3) = a_vec;
        dy.segment<4>(6) = pdot;
        dy.segment<4>(10) = pddot;
    }

    void constrain_rotations(state_type &y_in)
    {
        Map<VectorXd> y(y_in.data(), 14);
        y.segment<4>(6).normalize();
        double sigma = y.segment<4>(6).dot(y.segment<4>(10));
        y.segment<4>(10) -= sigma*y.segment<4>(6);
    }

    std::vector<VectorXd> integrate_to_target(
            double theta, // vertical angle in radians up is positive
            double phi, // horizontal angle in radians right is positive
            double x_target) // target x coordinate in meters
    {
        std::vector<VectorXd> result;

        double dy = tan(theta);
        double dz = tan(phi);
        Vector3d n_vec{1.0, dy, dz};
        n_vec.normalize();
        Vector3d z_vec = n_vec.cross(Vector3d{0.0, 1.0, 0.0});
        z_vec.normalize();
        Vector3d y_vec = z_vec.cross(n_vec);
        Matrix3d A_g;
        A_g.col(0) = n_vec;
        A_g.col(1) = y_vec;
        A_g.col(2) = z_vec;
        Vector4d p = caams::pA(A_g);
        double omega = k_rr*v0;
        Vector3d omega_p{omega, 0.0, 0.0};
        Vector4d pdot = caams::p_dot_omega_p(p, omega_p);

        state_type y(14);
        Map<VectorXd> y_map(y.data(), 14);
        y_map.segment<3>(0) = -y_vec*h_sight;
        y_map.segment<3>(3) = n_vec*v0;
        y_map.segment<4>(6) = p;
        y_map.segment<4>(10) = pdot;

        result.push_back(y_map);

        double t = 0.0;
        double delta_t = 2.0*M_PI/omega/32.0;
        double abs_err = 1e-16;
        double rel_err = 1e-16;
        auto stepper = make_controlled<error_stepper_type>(abs_err, rel_err);
        for(;;){
            integrate_adaptive(stepper, *this, y, 0.0, delta_t, delta_t);
            constrain_rotations(y);
            if(y[0] > x_target){
                break;
            }
            result.push_back(Map<VectorXd>(y.data(), 14));
        }

        double dx;
        while(fabs(dx = y[0]-x_target)>1e-4){
            double dt = dx/-y[3];
            integrate_adaptive(stepper, *this, y, 0.0, dt, dt);
            constrain_rotations(y);
        }

        result.push_back(Map<VectorXd>(y.data(), 14));

        return result;
    }

    void acquire_target(double &theta, double &phi, double x_target, double y_target)
    {
        for(;;){
            std::vector<VectorXd> result = integrate_to_target(theta, phi, x_target);
            VectorXd y_t = result.back();
            double dy = y_t(1) - y_target;
            double dz = y_t(2);
            if( fabs(dy)<1e-3 && fabs(dz)<1e-3 ){
                break;
            }
            double dtheta = atan(dy/x_target);
            double dphi = atan(dz/x_target);
            theta -= dtheta;
            phi -= dphi;
        }
    }
};

int main(void)
{
    Ballistics ballistics;
    double theta = 0.0;
    double phi = 0.0;
    double x_target = 100*m_p_yd;
    std::vector<VectorXd> result = ballistics.integrate_to_target(theta,phi,x_target);
    VectorXd y_t = result.back();
    std::cout.precision(17);
    std::cout << "Target coordinates:\n" << y_t.segment<3>(0) << std::endl;
    Vector4d p = y_t.segment<4>(6);
    Vector4d pdot = y_t.segment<4>(10);
    Vector3d omega_p = caams::omega_p_p_dot(p, pdot);
    double omega = omega_p.norm();
    std::cout << "omega target:" << omega << std::endl;
    Matrix3d A_g = caams::Ap(p);
    Vector3d x_vec = A_g.col(0);
    double moa_v = atan(x_vec(1)/x_vec(0))*moa_p_rad;
    double moa_h = atan(x_vec(2)/x_vec(0))*moa_p_rad;
    std::cout << "moa_v:" << moa_v << " moa_h:" << moa_h << std::endl;
    double theta_zero=0.0;
    double phi_zero=0.0;
    double x_zero = 200*m_p_yd;
    ballistics.acquire_target(theta_zero, phi_zero, x_zero, 0.0);
    std::cout << "theta_zero(MOA):" << theta_zero*moa_p_rad <<
                 " phi_zero(MOA):" << phi_zero*moa_p_rad << std::endl;
}



