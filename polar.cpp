/*
 * Integrates linear motion specified by polar
 * coordinates.
 *
 * state vector definition
 *
 * y[0] = theta
 * y[1] = dtheta
 * y[2] = radius
 * y[3] = dradius
 *
 * derivatives vector as a function of
 * state vector
 *
 * dy[0] = y[1]
 * dy[1] = -2*y[3]*y[1]/y[2]
 * dy[2] = y[3]
 * dy[3] = y[2]*y[1]*y[1]
 *
 */
#include <boost/numeric/odeint.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

using namespace boost::numeric::odeint;
using namespace Eigen;

//typedef VectorXd state_type;
typedef std::vector<double> state_type;


void dy_func(const state_type &y, state_type &dy, const double t)
{
    dy[0] = y[1];
    dy[1] = -2.0*y[3]*y[1]/y[2];
    dy[2] = y[3];
    dy[3] = y[2]*y[1]*y[1];
}


state_type initial_state(double x, double y, double vx, double vy)
{
    double radius = sqrt(x*x + y*y);
    double theta = atan2(y, x);
    double dr = vx*cos(theta) + vy*sin(theta);
    double dtheta = (vy*cos(theta) - vx*sin(theta))/radius;

    state_type r(4);
    r[0] = theta;
    r[1] = dtheta;
    r[2] = radius;
    r[3] = dr;
    return r;
}

void print_state(state_type y)
{
    double px = y[2]*cos(y[0]);
    double py = y[2]*sin(y[0]);
    std::cout.precision(17);
    std::cout << "<x,y>:"<< px << "," << py << std::endl;
}

int main(void)
{
    state_type y = initial_state(0.1, -10.0, 0.0, 10.0);

    typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
    controlled_stepper_type controlled_stepper;

    int N_points = 20;
    double t = 0.0;
    double dt = 2.0/N_points;
    double abs_err = 1e-16;
    double rel_err = 1e-16;
    print_state(y);
    for(int i=0;i<N_points;i++){
        size_t N_steps = integrate_adaptive(
                    make_controlled<error_stepper_type>(abs_err, rel_err),
                    dy_func , y , 0.0, dt , dt );
        std::cout << "N_steps:" << N_steps << std::endl;
        print_state(y);
        t += dt;
    }
    integrate_adaptive(
        make_controlled<error_stepper_type>(abs_err, rel_err),
        dy_func , y , 0.0, -dt , -dt );
    print_state(y);
    return 0;
}
