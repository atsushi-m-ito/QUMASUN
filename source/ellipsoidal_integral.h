
/********************************************************
* Ellipsoidal coordinate
*
* This is the 3 dimensional coordinate corresponding to Elliptic coordinate in 2 dimension.
*
* There are two foci F1 and F2 located at (0,0,-a) and (0,0,a) in the cartesian coordinate system (x,y,z).
* a = R/2.
* The r_1 and r_2 are the distance from the foci F1 and F2.
* The \phi is angle arround the z axis.
* Here, the Ellipsoidal coordinate (\xi,\eta,\phi) is given by
*   r_1 = (\xi + \eta)a
*   r_2 = (\xi - \eta)a
* then,
*   \xi = (r_1 + r_2)/2a
*   \eta = (r_1 - r_2)/2a
* The range of \xi and \eta are
*	\xi in [1, \infty]
*   \eta in [-1, 1]
* The Jacobian is as follows
*   dxdydz = a^3 (\xi^2 - \eta^2) d\xi d\eta d\phi
*
*********************************************************/

#include "GaussLegendre.h"
#include "GaussLaguerre.h"

namespace ReGZ {

    /*
    * Here, a and b is effective unit length that function roughly reduce as exp(-r_a/a) or exp(-r_b/b).
    * If you have no information of it, set a and b to 1.0.
    */
    template <class RET_TYPE, int N_xi, int N_eta, class FUNC>
    RET_TYPE IntegrateEllipsoidal(double R, double a, double b, FUNC func) {


        const double R_2 = R / 2.0; // R  is distance between two centers.
        const double dp = 1.0 / (1.0 / a + 1.0 / b);

        // range t in [0, \infty]
        RET_TYPE sum = Quadrature::Gauss_Laguerre::IntegrateGT<N_xi, RET_TYPE>(
            [&dp, &R_2, &func](double t) {
                const double xi = t * dp / R_2 + 1.0;

                const RET_TYPE sum_s = Quadrature::Gauss_Legendre::IntegrateT< N_eta, RET_TYPE>(
                    [&dp, &R_2, &xi, &func](double eta) {

                        const double r_a = (xi + eta) * R_2;
                        const double r_b = (xi - eta) * R_2;
                        const double J = r_a * r_b;
                        //const double J = xi * xi - eta * eta;

                        //const double cos_theta = (4.0 * R_2 * R_2 + r_a * r_a - r_b * r_b) / (4.0 * R_2 * r_a);
                        //const double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

                        return func(r_a, r_b) * J;
                    });
                return sum_s * dp;// / R_2;
            });

        return sum * (2.0 * M_PI);// *R_2* R_2* R_2;

    }

}
