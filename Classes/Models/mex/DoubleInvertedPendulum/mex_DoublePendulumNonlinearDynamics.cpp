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

#include "include/DoubleInvertedPendulum.h"
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include "armadillo"
#include "armaMex.hpp"
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;
using namespace boost::numeric::odeint;
using namespace Models;

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // mandatory arguments: sampling interval, massCart, massPendulum1, massPendulum2, lengthPendulum1, lengthPendulum2
    // frictionCart, frictionPendulum1, frictionPendulum2
    // state samples, input samples (might be empty), noise samples (might be empty)
    const double samplingInterval = armaGetDouble(inputArrays[0]);
    const double massCart = armaGetDouble(inputArrays[1]);
    const double massPendulum1 = armaGetDouble(inputArrays[2]);
    const double massPendulum2 = armaGetDouble(inputArrays[3]);
    const double lengthPendulum1 = armaGetDouble(inputArrays[4]);
    const double lengthPendulum2 = armaGetDouble(inputArrays[5]);
    const double frictionCart = armaGetDouble(inputArrays[6]);
    const double frictionPendulum1 = armaGetDouble(inputArrays[7]);
    const double frictionPendulum2 = armaGetDouble(inputArrays[8]);
    
    const dmat stateSamples = armaGetPr(inputArrays[9]);
    
    const uword numSamples = stateSamples.n_cols;
    // we return the states propagated in time
    outputArrays[0] = mxCreateUninitNumericMatrix(DoubleInvertedPendulum::DIM_STATE, static_cast<size_t>(numSamples), 
            mxDOUBLE_CLASS, mxREAL);
    dmat predictedStates(mxGetPr(outputArrays[0]), DoubleInvertedPendulum::DIM_STATE, numSamples, false, true);
    
    auto plant = new DoubleInvertedPendulum(massCart, massPendulum1, massPendulum2, lengthPendulum1,
                lengthPendulum2, frictionCart, frictionPendulum1, frictionPendulum2);
    
    const drowvec inputSamples = (!mxIsEmpty(inputArrays[10])) ? conv_to<drowvec>::from(armaGetPr(inputArrays[10])) : zeros<drowvec>(numSamples);
    const dmat noiseSamples = (!mxIsEmpty(inputArrays[11])) ? armaGetPr(inputArrays[11]) : zeros<dmat>(DoubleInvertedPendulum::DIM_NOISE, numSamples);
    
    #ifdef _OPENMP
    #pragma omp parallel for shared (predictedStates)
    #endif
    for (uword j=0; j < numSamples; ++j) {
        dvec6 x = stateSamples.unsafe_col(j);        
        const dvec3 w = noiseSamples.unsafe_col(j);
        plant->setInputForce(inputSamples(j));
        plant->setNoise(w);
        integrate_const(runge_kutta4<dcolvec6>(), *plant, x, 0.0, samplingInterval, samplingInterval);        
        predictedStates.col(j) = x;
    }
    delete plant;
}