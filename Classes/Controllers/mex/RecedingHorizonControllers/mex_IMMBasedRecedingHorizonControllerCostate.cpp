/* >> This file is part of CoCPN-Sim
*
*    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
*
*    Copyright (C) 2018-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // for computation of the costate, we need: terminalP, augA, augB, transitionMatrix, horizonLength, augR, augQ
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const dmat transitionMatrix = armaGetPr(inputArrays[2]);
    const dcube augQ = armaGetCubePr(inputArrays[3]);
    const dmat augR = armaGetPr(inputArrays[4]); // only for first mode
    const dcube terminalP = armaGetCubePr(inputArrays[5]);   
    const int horizonLength = armaGetScalar<int>(inputArrays[6]);
    
    const uword numModes = augA.n_slices;    
    const uword dimAugState = augA.n_rows;
     
    size_t dims[3] = {static_cast<size_t>(dimAugState), static_cast<size_t>(dimAugState), static_cast<size_t>(numModes)};
    outputArrays[0] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // the costate    
    dcube P(mxGetPr(outputArrays[0]), dimAugState, dimAugState, numModes, false, true);    
    P = terminalP;
    const dmat eyeM = eye<dmat>(augB.n_cols, augB.n_cols);
    
    for (auto k = horizonLength - 1; k >=1; --k) {
        // backward iteration over the horizon
        const dcube P_old = P;
        
        #ifdef _OPENMP        
        #pragma omp parallel for shared(P)        
        #endif
        for (uword i = 0; i < numModes; ++i) {
            // compute Y_k+1 using P_k+1
            dmat Y = zeros<dmat>(dimAugState, dimAugState);
            mat::const_row_iterator row = transitionMatrix.begin_row(i);
            for (uword j = 0; j < numModes; ++j, ++row) {
                Y += (*row) * P_old.slice(j);
            }                
        
            // fill M
            const dmat BY = augB.slice(i).t() * Y;
            // this is the expression for M in the paper                 
            const dmat M = (i== 0) ? conv_to<dmat>::from(BY * augB.slice(i) + augR) : conv_to<dmat>::from(BY * augB.slice(i)); // symmetric
            const dmat YBMBY = symmatu(Y - BY.t() * solve(M, eyeM) * BY);                        
            P.slice(i) = symmatu(augQ.slice(i) + augA.slice(i).t() * YBMBY * augA.slice(i)); // symmetric            
        }           
    }    
}
