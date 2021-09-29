/* >> This file is part of CoCPN-Sim
*
*    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
*
*    Copyright (C) 2018-2020  Florian Rosenthal <florian.rosenthal@kit.edu>
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
 
#include "armadillo"
#include "armaMex.hpp"
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

const double GRAVITATIONAL_ACCELERATION = 9.81;
const int NUM_STATES = 4; // cart position, velocity, pendulum angle, angular velocity

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // mandatory arguments: sampling interval, inertia, massPendulum, massCart, lengthPendulum, friction
    // state samples, input samples (might be empty), noise samples (might be empty)
    const double samplingInterval = armaGetDouble(inputArrays[0]);
    const double inertia = armaGetDouble(inputArrays[1]);
    const double massPendulum = armaGetDouble(inputArrays[2]);
    const double massCart = armaGetDouble(inputArrays[3]);
    const double lengthPendulum = armaGetDouble(inputArrays[4]);
    const double friction = armaGetDouble(inputArrays[5]);
    const dmat stateSamples = armaGetPr(inputArrays[6]);
        
    const uword numSamples = stateSamples.n_cols;
    // we return the states propagated in time
    outputArrays[0] = mxCreateUninitNumericMatrix(NUM_STATES, static_cast<size_t>(numSamples), 
            mxDOUBLE_CLASS, mxREAL);
    
    dmat predictedStates(mxGetPr(outputArrays[0]), NUM_STATES, numSamples, false, true);
    
    // new cart position (x_1) and pendulum angle (x_3)   
    predictedStates.row(0) = samplingInterval * stateSamples.row(1) + stateSamples.row(0);
    predictedStates.row(2) = samplingInterval * stateSamples.row(3) + stateSamples.row(2);
    
    // compute new cart velocity (x_2) and pendulum angular velocity (x_4)
    const double massSum = massCart + massPendulum; // m + M
    const double massLength = massPendulum * lengthPendulum; // m * l
    const double inertiaSum = inertia + massLength * lengthPendulum; // I+m*l^2
    const drowvec sines = sin(stateSamples.row(2));
    //const drowvec cosines = cos(stateSamples.row(2));
    const drowvec massLengthCos = massLength * cos(stateSamples.row(2)); // m*l* cos(x_3)
    const drowvec denominator = massSum * inertiaSum - square(massLengthCos); // (M+m)*(I+m*l^2)-m^2*l^2*cos(x_3)^2
    // u + m*l*sin(x_3)*(x_4)^2-b*x_2 + w_u
    // noisy input here (actuation forces enters the same way as control input)
    drowvec forces = massLength * sines % square(stateSamples.row(3)) - friction * stateSamples.row(1);
    if (!mxIsEmpty(inputArrays[7])) {
        // add the input samples
        forces += conv_to<drowvec>::from(armaGetPr(inputArrays[7]));
    }
    // additional forces: m * g * sin(x_3) - w_p
    // disturbance force here (disturbance force acting on the pendulum tip)
    drowvec forces2 =  massPendulum * GRAVITATIONAL_ACCELERATION * sines;
    if (!mxIsEmpty(inputArrays[8])) {
        // add disturbance forces w_p and w_u
        const dmat disturbanceForces = armaGetPr(inputArrays[8]);
        forces2 -= disturbanceForces.row(0);
        forces += disturbanceForces.row(1);
    }
    predictedStates.row(1) = stateSamples.row(1) 
        + samplingInterval * (inertiaSum * forces + lengthPendulum * massLengthCos % forces2) / denominator;
    predictedStates.row(3) = stateSamples.row(3) 
        + samplingInterval * (massSum * lengthPendulum * (-forces2) - massLengthCos % forces) / denominator;  
}
