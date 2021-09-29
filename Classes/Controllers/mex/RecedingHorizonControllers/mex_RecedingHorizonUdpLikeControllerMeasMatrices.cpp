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

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {    
    // mandatory arguments: modesScBase2, C, dimEta
    const umat modesScBase2 = conv_to<umat>::from(armaGetPr(inputArrays[0]));
    const dmat C = armaGetPr(inputArrays[1]); // base meas matrix C
    const uword dimEta = armaGetScalar<uword>(inputArrays[2]);
    
    const uword measBufferLength =  modesScBase2.n_cols;
    const uword numModesSc = modesScBase2.n_rows;    
    const uword dimMeas = C.n_rows;
    const uword dimAugmentedMeas = dimMeas * measBufferLength;
    
    size_t dims[3] = {static_cast<size_t>(dimAugmentedMeas), static_cast<size_t>(dimAugmentedMeas), static_cast<size_t>(numModesSc)};
    outputArrays[0] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // S matrices   
    dcube S(mxGetPr(outputArrays[0]), dimAugmentedMeas, dimAugmentedMeas, numModesSc, false, true);
    
    const dmat mOnes = ones<drowvec>(dimMeas);
    #ifdef _OPENMP    
    #pragma omp parallel for shared (S)
    #endif
    for (uword j=0; j < numModesSc; ++j) {
        S.slice(j) = diagmat(kron(conv_to<dmat>::from(modesScBase2.row(j) == 1), mOnes)); 
    }
    
    size_t dimsAugC[2] = {static_cast<size_t>(dimAugmentedMeas), static_cast<size_t>(measBufferLength * C.n_cols + dimEta)};
    sp_dmat augC(dimsAugC[0], dimsAugC[1]);   
   
    // populate
    for (uword i=0; i < measBufferLength; ++i) {
         augC.submat(i * dimMeas, i * C.n_cols, (i + 1) * dimMeas - 1, (i + 1) * C.n_cols - 1) = C; 
    }    
    //augC.head_cols(measBufferLength * C.n_cols) = kron(eye(measBufferLength, measBufferLength), C);    
        
    outputArrays[1] = mxCreateSparse(dimsAugC[0], dimsAugC[1], augC.n_nonzero, mxREAL); // augC , sparse   
    armaSetSparsePr(outputArrays[1], augC);   
}

