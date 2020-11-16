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
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {    
    // mandatory arguments: modesScBase3, measBufferLength, normalizedProbs 
    const umat modesScBase3 = conv_to<umat>::from(armaGetPr(inputArrays[0]));
    const uword measBufferLength = armaGetScalar<uword>(inputArrays[1]);
    const drowvec normalizedProbs = armaGetPr(inputArrays[2]);
    
    const uword numModesSc = modesScBase3.n_rows;   
    const uword maxMeasDelay = measBufferLength - 1;
    const drowvec sums = cumsum(normalizedProbs);
    const dmat33 trans = {  {0, 0, 0},
                            {0, 0, 1},
                            {0, 0, 1},
                         };
    
    dcube measStateTransitionMatrices(3, 3, measBufferLength);
    measStateTransitionMatrices.each_slice() = trans;
    
    const dvec a_i = normalizedProbs.subvec(1, maxMeasDelay + 1) / (1 - sums.head(measBufferLength));
    measStateTransitionMatrices.tube(0,0) = 1 - a_i;
    measStateTransitionMatrices.tube(0,1) = a_i;
    
    const uword numRows = pow(3, maxMeasDelay);
    outputArrays[0] = mxCreateNumericMatrix(static_cast<size_t>(numRows), static_cast<size_t>(numModesSc), mxDOUBLE_CLASS, mxREAL);    
    dmat transitionMatrixSc_part(mxGetPr(outputArrays[0]), numRows, numModesSc, false, true);
    
    #ifdef _OPENMP    
    #pragma omp parallel for shared (transitionMatrixSc_part)
    #endif
    for (uword j = 0; j < numModesSc; ++j) {
        const uvec j_idx = modesScBase3.row(j); // new mode
        double* col = transitionMatrixSc_part.colptr(j);
        for (uword i = 0; i < numRows; ++i) {
            const uvec i_idx = modesScBase3.row(i); // old mode
            double prob = (j_idx.at(0) == 0) ? (1 - normalizedProbs.at(0)) : ((j_idx.at(0) == 1) ? normalizedProbs.at(0) : 0);
            uword a = 0;
            while (prob != 0 && a < maxMeasDelay) {
                prob *= measStateTransitionMatrices.at(i_idx(a), j_idx(a+1), a);
                ++a;
            }
            col[i] = prob;
        }
    }
}