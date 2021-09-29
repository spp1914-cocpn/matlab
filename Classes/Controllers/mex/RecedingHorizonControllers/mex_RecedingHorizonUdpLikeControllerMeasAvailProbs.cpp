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
 
#include "armadillo"
#include "armaMex.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

double evalDepth(const uword depth, const uword state, const uword colIdx, const field<uvec>& measAvailabilityIdx, 
        const dmat& scDelayTransMat, const dcolvec& delayProbs) {
    uvec localIdx = measAvailabilityIdx(depth, state);    
    if (depth == measAvailabilityIdx.n_rows -1) {            
        // portion of T
        const vec rows = scDelayTransMat.unsafe_col(colIdx);
        
        return accu(rows.elem(localIdx) % delayProbs.elem(localIdx));        
        // final stage of recursion: T(m,j) * p(m) (last element of continued expression)
    }
    
    double prob = 0;
    for (auto rowIdx : localIdx) {
        prob += scDelayTransMat(rowIdx-1, colIdx) 
            * evalDepth(depth + 1, state, rowIdx, measAvailabilityIdx, scDelayTransMat, delayProbs);
        // enter recursion: T(i, s) * [T(l, i) * .....T(j, r) * T(m,j) * p(m)]
    }
    return prob;
}

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {    
    // mandatory arguments: sc delay transition matrix, meas availability indices (cell array!), delay probs (column vector)
    const dmat scDelayTransMat = armaGetPr(inputArrays[0]);
    const dcolvec delayProbs = armaGetPr(inputArrays[2]);  
    
    // measAvailabilityIdx = cell(measBufferLength, 2^measBufferLength), one column per state
    // move content into a field    
    const mwSize* cellDims = mxGetDimensions(inputArrays[1]);
    const uword measBufferLength = static_cast<uword>(cellDims[0]); // rows
    const uword numStates = static_cast<uword>(cellDims[1]); // columns
    field<uvec> measAvailabilityIdx(measBufferLength, numStates);
    
    mwIndex matIdx = 0;
    for (uword i=0; i < numStates; ++i) {
        for (uword j=0; j < measBufferLength; ++j) {
           uvec indices = conv_to<uvec>::from(armaGetPr(mxGetCell(inputArrays[1], matIdx++))); // indices start at 1           
           indices -= 1; //0 based indices in C/C++
           measAvailabilityIdx(j, i) = indices;
        }
    }
 
    outputArrays[0] = mxCreateDoubleMatrix(1, static_cast<size_t>(numStates), mxREAL); // row vector, filled with zeros    
    double* probs = static_cast<double*>(mxGetPr(outputArrays[0]));
    
    switch (measBufferLength) {
        case 1:
            // corner case, no delayed measurement allowed
            probs[0] = accu(delayProbs.elem(measAvailabilityIdx(0, 0)));
            probs[1] = accu(delayProbs.elem(measAvailabilityIdx(0, 1)));
            break;
        case 2:
            // max meas delay is 1
            for (auto k=0; k < 4; ++k) {
                const uvec localIdx = measAvailabilityIdx(1, k);                
                for (auto l : measAvailabilityIdx(0, k)) {
                    probs[k] += accu(scDelayTransMat.unsafe_col(l).elem(localIdx) % delayProbs.elem(localIdx));
                }
            }
            break;
        case 3:
            // max meas delay is 2
            for (auto k=0; k < 8; ++k) {
                const uvec localIdx = measAvailabilityIdx(2, k);
                const uvec localIdx2 = measAvailabilityIdx(1, k);
                for (auto l : measAvailabilityIdx(0, k)) {
                    for (auto m : localIdx2) {
                        probs[k] += scDelayTransMat(m, l) * accu(scDelayTransMat.unsafe_col(m).elem(localIdx) % delayProbs.elem(localIdx));
                    }
                }                
            }
            break;
        default:
            // compute probabilities recursively
            for (uword k=0; k < numStates; ++k) {
                const uvec currIdx = measAvailabilityIdx(0, k);
                //double* probs = static_cast<double*>(mxGetPr(outputArrays[0]));
                for (auto l : currIdx) {            
                    probs[k] += evalDepth(1, k, l, measAvailabilityIdx, scDelayTransMat, delayProbs);             
                }        
            }
            break;
    }   
}
