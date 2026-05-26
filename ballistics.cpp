#include "conversions.h"
#include "caams.hpp"
#include <boost/numeric/odeint.hpp>
// #include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
// #include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
// #include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
// #include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_fehlberg78.hpp>
// #include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <mutex>
#include <thread>
#include <chrono>

using namespace std::chrono_literals;
using namespace boost::numeric::odeint;
using namespace Eigen;
typedef std::vector<double> state_type;
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;

class Ballistics
{
public:
    // Global
    double g;
    double omega_earth;
    Vector3d omega_earth_vec;

    // Projectile
    double m;
    double v0;
    double diameter;
    double length;
    double d_cp;
    double radius;
    double A;
    Matrix3d Jp;
    Vector3d r_cp;
    double theta_zero; // vertical zero offset
    double phi_zero;   // horizontal zero offset
    double range_zero; // range for zero settings

    // Rifle
    double d_rr;
    double k_rr;
    double h_sight;

    // Range
    double heading_range;
    double altitude;
    double latitude;
    Vector3d r_base;
    Matrix3d A_nue; // north up east
    Matrix3d A_base; // range orientation
    Matrix3d A_target; // target orientation
    Vector3d r_target;
    Vector3d N_target; // target plane

    // Atmosphere
    double heading_wind;
    double vw;
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

    void r_base_init(void)
    {
        const double e2 = 0.0066943800229;
        const double a = 6378137.0;
        double phi = latitude*M_PI/180.0;
        double sin_phi = sin(phi);
        double cos_phi = cos(phi);
        double N = a/sqrt(1-e2*sin_phi*sin_phi);
        r_base << (N+altitude)*cos_phi,
                  0.0,
                (N*(1-e2)+altitude)*sin_phi;
    }

    void A_init(void)
    {
        r_base_init();
        Vector3d y_axis = r_base.normalized();
        Vector3d east = Vector3d{0.0, 0.0, 1.0}.cross(y_axis).normalized();
        Vector3d north = y_axis.cross(east);
        A_nue.col(0) = north;
        A_nue.col(1) = y_axis;
        A_nue.col(2) = east;
        double theta = heading_range*M_PI/180;
        Vector3d x_axis = cos(theta)*north +
                sin(theta)*east;
        Vector3d z_axis = x_axis.cross(y_axis);
        A_base.col(0) = x_axis;
        A_base.col(1) = y_axis;
        A_base.col(2) = z_axis;
    }

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

    void wind_vec_init(void)
    {
        A_init();
        double theta = heading_wind*M_PI/180.0;
        vw_vec = A_nue * Vector3d{
                cos(theta),
                0.0,
                sin(theta)} * vw;
    }

    double Cd_at_mach(double m)
    {
        for(int i=0;;i++){
            if( m>=Cd_table(0,i) && m<=Cd_table(0,i+1)){
                double alpha = (m-Cd_table(0,i))/(Cd_table(0,i+1)-Cd_table(0,i));
                return Cd_table(1,i)*(1.0-alpha) + Cd_table(1,i+1)*alpha;
            }
        }
    }

    double speed_of_sound(double T, double p, double humid)
    {
        double pv = humid/100.0*psat(T);
        double pd = p - pv;
        double rho = rho_from_weather(T,p,humid);
        double K = 1.3*pv + 1.4*pd;
        return sqrt(K/rho);
    }

    Ballistics(void):
        g(9.80665),
        omega_earth(2.0*M_PI/86164.0905),
        omega_earth_vec(0.0, 0.0, omega_earth),
        // Projectile
        m(120.0*kg_p_gr),
        v0(2900.0*mps_p_fps),
        diameter(6.7e-3),
        length(60.0e-3),
        d_cp(10.0e-3),
        radius(diameter/2),
        A(M_PI*radius*radius),
        h_sight(1.5*m_p_inch),
        Jp(caams::J_p_cylinder_x_axis(m, radius, length)),
        r_cp(d_cp, 0.0, 0.0),

        // barrel
        d_rr(8*m_p_inch),
        k_rr(2*M_PI/d_rr),

        // atmospheric parameters
        vw(0*mps_p_kmph),
        heading_range(90.0),
        altitude(1000.0),
        latitude(0.0),
        heading_wind(0.0),

        R(8.31466),
        Md(0.0289652),
        Mv(0.018016),
        L(0.0065),
        p0(29.53*Pa_p_inHg),
        Tc(0.0),
        T0(Tc + 273.15),
        humidity(78),
        rho(rho_from_weather(T0, p0, humidity)),
        c(speed_of_sound(T0, p0, humidity))
    {
        Cd_table.resize(2,9);
        Cd_table << 0.,   0.7,  0.8,  0.9,  1.0, 1.1,   1.6,  3.5,  5.0,
                    0.12, 0.12, 0.13, 0.15, 0.4, 0.405, 0.33, 0.22, 0.17;

        wind_vec_init();

        printf("Speed of sound:%.1lf\n", c);
    }

    //
    // This class is a functor and the function computes the derivatives
    // of the state vector.
    //
    void operator() (const state_type &y_in, state_type &dy_in, const double t)
    {
        Map<const VectorXd> y(y_in.data(), 14);
        Map<VectorXd> dy(dy_in.data(), 14);
        Vector3d r = y.segment<3>(0);
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
        Vector3d a_g_vec = (-r).normalized()*g;
        // coriolis force
        Vector3d a_coriolis = -2.0*omega_earth_vec.cross(v_vec);
        // centrifugal force
        Vector3d a_cent = -omega_earth_vec.cross(
                    omega_earth_vec.cross(r));
        Vector3d a_vec = a_g_vec + a_coriolis + a_cent + F_air/m;
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

    std::vector<VectorXd> result;
    std::mutex result_mtx;
    double d_target;
    bool thread_finished;

    void result_clear(void){
        std::lock_guard<std::mutex> lock(result_mtx);
        result.clear();
    }

    void result_push(state_type &y)
    {
        std::lock_guard<std::mutex> lock(result_mtx);
        result.push_back(Map<VectorXd>(y.data(), 14));
    }

    VectorXd result_last(void)
    {
        //printf("result_last. testing lock\n");
        std::lock_guard<std::mutex> lock(result_mtx);
        //printf("mutex acquired.\n");
        if(result.empty()){
            return VectorXd::Zero(14);
        }else{
            return result.back();
        }
    }

    std::vector<VectorXd> result_all(void)
    {
        std::lock_guard<std::mutex> lock(result_mtx);
        std::vector<VectorXd> new_result = result;
        return new_result;
    }

    double percent_done(void){
        VectorXd y = result_last();
        Vector3d r_g = y.segment<3>(0);
        Vector3d r_p = r_g - r_base;
        double d = r_p.dot(A_target.col(0));
        return d/d_target*100.0;
    }

    void integrate_to_target(
            double theta,
            double phi,
            double target_elevation,
            double target_distance)
    {
        result_clear();
        thread_finished = false;
        d_target = target_distance;

        printf("integrate_to_target\n");

        std::thread thd(
                    &Ballistics::integrate_to_target_thread_func, this,
                    theta, phi, target_elevation, target_distance);

        while(!thread_finished){
            double p = percent_done();
            printf("Percent done: % 6.1lf%%\r", p);
            fflush(stdout);
            std::this_thread::sleep_for(25ms);
        }
        printf("\n");
        thd.join();
    }

    void integrate_to_target_thread_func(
            double theta, // vertical angle in radians up is positive
            double phi, // horizontal angle in radians right is positive
            double target_elevation,
            double target_distance) //
    {
        A_init();
        // Target basis. Rotate A_base about A_base.col(2)
        double angle = target_elevation*M_PI/180.0;
        Matrix3d A_rot = caams::AAA(angle, A_base.col(2));
        A_target = A_rot*A_base;
        // target location and normal
        r_target = r_base + A_target.col(0)*target_distance;
        N_target = -A_target.col(0);
        // rifle (projectile) basis including aim offsets
        double dy = tan(theta);
        double dz = tan(phi);
        Vector3d n_vec{1.0, dy, dz};
        Vector3d x_rifle = A_target*n_vec;
        x_rifle.normalize();
        Vector3d z_rifle = x_rifle.cross(A_target.col(1)).normalized();
        Vector3d y_rifle = z_rifle.cross(x_rifle);
        Matrix3d A_rifle;
        A_rifle.col(0) = x_rifle;
        A_rifle.col(1) = y_rifle;
        A_rifle.col(2) = z_rifle;
        // projectile orientation quaternion and derivative
        Vector4d p = caams::pA(A_rifle);
        double omega = k_rr*v0;
        Vector3d omega_p{omega, 0.0, 0.0};
        Vector4d pdot = caams::p_dot_omega_p(p, omega_p);

        state_type y(14);
        Map<VectorXd> y_map(y.data(), 14);
        y_map.segment<3>(0) = r_base-h_sight*y_rifle;
        y_map.segment<3>(3) = x_rifle*v0;
        y_map.segment<4>(6) = p;
        y_map.segment<4>(10) = pdot;

        double t = 0.0;
        double delta_t = 2.0*M_PI/omega/32.0;
        double delta_t_step;
        double abs_err = 1e-16;
        double rel_err = 1e-16;
        auto stepper = make_controlled<error_stepper_type>(abs_err, rel_err);
        for(;;){
            result_push(y);
            // estimate impact time
            Map<VectorXd> ym(y.data(), 14);
            Vector3d r0 = ym.segment<3>(0);
            Vector3d v = ym.segment<3>(3);
            Vector3d delta_r = r0 - r_target;
            double dist = delta_r.dot(N_target);
            if(fabs(dist)<1e-3){
                break;
            }
            double t_impact = -dist/
                    v.dot(N_target);
            // select time step
            if(fabs(t_impact)<delta_t){
                delta_t_step = t_impact;
            }else{
                delta_t_step = delta_t;
            }
            integrate_adaptive(stepper,
                               std::ref(*this), y, 0.0,
                               delta_t_step, delta_t_step);
            constrain_rotations(y);
        }

        thread_finished = true;
    }

    void acquire_target(double &theta, double &phi,
                        double target_elevation,
                        double target_distance)
    {
        for(;;){
            integrate_to_target(theta, phi,
                                target_elevation,
                                target_distance);
            VectorXd y_t = result_last();
            Vector3d r_hit_g = y_t.segment<3>(0) - r_base;
            // resolve the hit in the target space
            Vector3d r_t = A_target.transpose()*r_hit_g;
            double dy = r_t(1);
            double dz = r_t(2);
            if( fabs(dy)<1e-3 && fabs(dz)<1e-3 ){
                break;
            }
            double dtheta = atan(dy/target_distance);
            double dphi = atan(dz/target_distance);
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
    double target_distance = 200*m_p_yd;
    double target_elevation = 0.0;
    ballistics.integrate_to_target(
                theta,phi,
                target_elevation,
                target_distance);
    VectorXd y_t = ballistics.result_last();
    Vector3d r_hit_g = y_t.segment<3>(0) - ballistics.r_base;
    // resolve the hit in the target space
    Vector3d r_t = ballistics.A_target.transpose()*r_hit_g;
    double dy = r_t(1);
    double dz = r_t(2);

    std::cout.precision(17);
    std::cout << "Target coordinates:\n" <<
                 "vertical(inch):" << dy/m_p_inch <<
                 " spin drift(inch):" << dz/m_p_inch << std::endl;

    target_distance = 200.0*m_p_yd;
    target_elevation = 0.0;
    ballistics.acquire_target(theta, phi, target_elevation,
                              target_distance);

}



