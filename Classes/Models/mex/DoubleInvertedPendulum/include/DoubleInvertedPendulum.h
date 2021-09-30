/* >> This file is part of CoCPN-Sim
*
*    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
*
*    Copyright (C) 2018-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
*
*                        Institute for Anthropomatics and Robotics
*                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
*                        Karlsruhe Institute of Technology (KIT), Germany
*
*                        https://isas.iar.kit.edu
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _MODELS_DOUBLE_INVERTED_PENDULUM_H
#define _MODELS_DOUBLE_INVERTED_PENDULUM_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <armadillo>
#include <cmath>

using namespace arma;

namespace Models { 
    
    class DoubleInvertedPendulum {
        
        private:
            double massCart; // positive scalar denoting the mass (in kg) of the cart
            double massPendulum1; // positive scalar denoting the mass of the bob (in kg) located at the top of the lower pendulum
            double massPendulum2; // positive scalar denoting the mass of the bob (in kg) located at the top of the upper pendulum
            double lengthPendulum1; // positive scalar denoting the length the lower pendulum rod (in m)
            double lengthPendulum2; // positive scalar denoting the length the upper pendulum rod (in m)
            double frictionCart; // positive scalar denoting the coefficient of viscous friction for the cart (in Ns/m)
            double frictionPendulum1; // nonnegative scalar denoting the coefficient of viscous friction for the joint at the lower pendulum (in N*s*m)
            double frictionPendulum2; // nonnegative scalar denoting the coefficient of viscous friction for the joint at the upper pendulum (in N*s*m).
            
            double inputForce; // input u is the force (in N) applied to move the cart, default 0
            dvec3 noise; //  w=[w1 w2 w3] is composed of disturbances (torques) acting at the tips of both pendulum rods (w1 and w2) and a disturbance force (w3) acting on the cart (additive on u, actuation noise), default [0; 0; 0]
            
            static constexpr double GRAVITATIONAL_ACCELERATION = 9.81;
        
        public:
            static constexpr double DIM_STATE = 6; // cart position, pendulum angle 1, pendulum angle 2 and corresponding velocities
            static constexpr double DIM_NOISE = 3; // noise entering as disturbance force at the cart, disturbances (torques) at either joint
            
        public:
            DoubleInvertedPendulum(const double massCart, const double massPendulum1, const double massPendulum2, const double lengthPendulum1,
                    const double lengthPendulum2, const double frictionCart, const double frictionPendulum1, const double frictionPendulum2,
                    const double inputForce = 0, const dvec3& noise = zeros<dcolvec>(3)) 
                : massCart(massCart), massPendulum1(massPendulum1), massPendulum2(massPendulum2), 
                    lengthPendulum1(lengthPendulum1), lengthPendulum2(lengthPendulum2), frictionCart(frictionCart), 
                    frictionPendulum1(frictionPendulum1), frictionPendulum2(frictionPendulum2),
                    inputForce(inputForce), noise(noise) {}   
            
            DoubleInvertedPendulum(const DoubleInvertedPendulum& other) 
                : massCart(other.massCart), massPendulum1(other.massPendulum1), massPendulum2(other.massPendulum2), 
                    lengthPendulum1(other.lengthPendulum1), lengthPendulum2(other.lengthPendulum2), frictionCart(other.frictionCart), 
                    frictionPendulum1(other.frictionPendulum1), frictionPendulum2(other.frictionPendulum2),
                    inputForce(other.inputForce), noise(other.noise) {}  
  
            
            ~DoubleInvertedPendulum() {}
            
            DoubleInvertedPendulum& operator=(const DoubleInvertedPendulum& other) {
                if (this == &other) {
                    return *this;
                }    
                this->massCart = other.massCart;
                this->massPendulum1 = other.massPendulum1;
                this->massPendulum2 = other.massPendulum2;
                this->lengthPendulum1 = other.lengthPendulum1;
                this->lengthPendulum2 = other.lengthPendulum2;
                this->frictionCart = other.frictionCart;           
                this->frictionPendulum1 = other.frictionPendulum1;
                this->frictionPendulum2 = other.frictionPendulum2;
                this->inputForce = other.inputForce;
                this->noise = other.noise;
                
                return *this;
            }  
            
            DoubleInvertedPendulum& operator=(DoubleInvertedPendulum&& other) {
                if (this == &other) {
                    return *this;
                }
                this->massCart = std::exchange(other.massCart, 0);
                this->massPendulum1 = std::exchange(other.massPendulum1, 0);
                this->massPendulum2 = std::exchange(other.massPendulum2, 0);
                this->lengthPendulum1 = std::exchange(other.lengthPendulum1, 0);
                this->lengthPendulum2 = std::exchange(other.lengthPendulum2, 0);     
                this->frictionCart = std::exchange(other.frictionCart, 0);           
                this->frictionPendulum1 = std::exchange(other.frictionPendulum1, 0);
                this->frictionPendulum2 = std::exchange(other.frictionPendulum2, 0);
                this->inputForce = std::exchange(other.inputForce, 0);
                this->noise = std::move(other.noise);

                return *this;
         }
        
        void operator() (const dcolvec6& x, dcolvec6& dxdt, const double t) {
            /*
             * dynamics is of the form M(q)q_dotdot = f(q,q_dot, u, w) (1)
             * with q the generalized coordinates cart position (in m),
             * angular deviation from upright for both pendulums (in rad)
             * q_dot and q_dotdot are the corresponding time derivatives
             * (i.e., generalized velocities and accelerations)
             * u is the input force (in N) applied to move the cart
             * w are the disturbance forces acting on the cart and both pendulum rods
             * M is the nonsingular mass matrix for all possible configurations q
             * f contains coriolis and gravity terms
             *
             * (1) can be rewritten in ODE form:
             * [I 0; 0 M]x_dot = [q_dot; f(x, u, w)] (2)
             * by introducing the state variable x = [q, q_dot]', hence x_dot = [q_dot, q_dotdot]
             * (2) can be written as x_dot = [q_dot; inv(M)*f(x,u,w)]
             * assuming that u and w are constant over the sampling interval, we can solve this ODE numerically
            */             
            const double massPend = this->massPendulum1 + this->massPendulum2;
            dmat invM = this->computeInverseMassMatrix(x.head(3)); // first three elements of x are q 
            
            const double sinDiff = sin(x(1) - x(2));
                        
            dxdt.head(3) = x.tail(3); // q_dot = q_dot
            
            // coriolis and friction terms
            dxdt(3) = -this->frictionCart * x(3) + this->lengthPendulum1 * massPend * sin(x(1)) * std::pow(x(4), 2) 
                + this->massPendulum2 * this->lengthPendulum2 * sin(x(2)) * std::pow(x(5), 2);
            dxdt(4) = -this->frictionPendulum1 * x(4) - this->lengthPendulum1 * this->lengthPendulum2 
                * this->massPendulum2 * sinDiff * std::pow(x(5), 2);
            dxdt(5) = -this->frictionPendulum2 * x(5) + this->lengthPendulum1 * this->lengthPendulum2 * this->massPendulum2 
                    * sinDiff * std::pow(x(4), 2);
            
            // gravity terms
            dxdt(4) += DoubleInvertedPendulum::GRAVITATIONAL_ACCELERATION * this->lengthPendulum1 * massPend * sin(x(1));
            dxdt(5) += DoubleInvertedPendulum::GRAVITATIONAL_ACCELERATION * this->lengthPendulum2 * this->massPendulum2 * sin(x(2));
            
            // noise and input
            dxdt(3) += as_scalar(this->inputForce + this->noise.tail(1));
            dxdt.tail(2) += this->noise.head(2);
            
            // lhs now contains q_dotdot
            dxdt.tail(3) = invM * dxdt.tail(3);
        }
    
        void setInputForce(const double inputForce) {
            // input force u, assumed constant over sampling interval
            this->inputForce = inputForce;
        }
        
        void setNoise(const dvec3& noise) {
            // disturbance/noise w, assumed constant over sampling interval
            this->noise = noise;
        }
        
    private:
        dmat33 computeInverseMassMatrix(const dvec3& q) {
            // symmetric mass matrix M and its inverse only dependent on pendulum angles
            const double massPend = this->massPendulum1 + this->massPendulum2;
            const double angleDiff = q(1)-q(2);
            
            const double detM = std::pow(this->lengthPendulum1 * this->lengthPendulum2, 2) * this->massPendulum2
                * (this->massCart * this->massPendulum1 + std::pow(this->massPendulum1 * sin(q(1)), 2)
                + this->massCart * this->massPendulum2 * std::pow(sin(angleDiff), 2)
                + this->massPendulum1 * this->massPendulum2 * std::pow(sin(q(1)), 2));
            
            // entries of the symmetric mass matrix M (state-dependent)
            const double M11 = this->massCart + massPend;
            const double M12 = this->lengthPendulum1 * massPend * cos(q(1));
            const double M13 = this->massPendulum2 * this->lengthPendulum2 * cos(q(2));
            const double M22 = std::pow(this->lengthPendulum1, 2) * massPend;
            const double M23 = this->lengthPendulum1 * this->lengthPendulum2 * this->massPendulum2 * cos(angleDiff);
            const double M33 = std::pow(this->lengthPendulum2, 2) * this->massPendulum2;
            
            return {    {(M22 * M33 - std::pow(M23, 2)) / detM , (M13 * M23 - M12 * M33) / detM, (M12 * M23 - M13 * M22) / detM},
                        {(M13 * M23 - M12 * M33) / detM, (M11 * M33 - std::pow(M13, 2)) / detM, (M13 * M12 - M11 * M23)/ detM},
                        {(M12 * M23 - M13 * M22) / detM, (M13 * M12 - M11 * M23) / detM, (M11 * M22 - std::pow(M12, 2)) / detM} };
            
        }
    };
}
#endif
