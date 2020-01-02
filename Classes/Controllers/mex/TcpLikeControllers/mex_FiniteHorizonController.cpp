/* >> This file is part of CoCPN-Sim
*
*    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
*
*    Copyright (C) 2017-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
#include <armaMex.hpp>

using namespace std;
using namespace arma;

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // required params: augA, augB, augQ, augR, transitionMatrix, terminalK, horizonLength
            
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const dcube augQ = armaGetCubePr(inputArrays[2]);
    const dcube augR = armaGetCubePr(inputArrays[3]);
    const dmat transitionMatrix = armaGetPr(inputArrays[4]);
    const dmat terminalK = armaGetPr(inputArrays[5]);
    const int horizonLength = armaGetScalar<int>(inputArrays[6]);
    
    const uword numModes = augA.n_slices;
    
    size_t dims[4] = {augR.n_rows, augA.n_rows, numModes, static_cast<size_t>(horizonLength)};
    outputArrays[0] = mxCreateUninitNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL); // the gains L
    if (outputArrays[0] == NULL) {
        // something went wrong
        mexErrMsgIdAndTxt("FiniteHorizonController:mex_FiniteHorizonController", 
                "** Allocation of output array failed **");
        return;
    }
    
    // we only need K_{k+1} to compute K_k and L_k
    dcube K_k(size(augA)); 
    K_k.each_slice() = terminalK; // K_N is the same for all modes
        
    double* dst = mxGetPr(outputArrays[0]);
    dst += ((horizonLength-1) * augR.n_rows * augA.n_rows * numModes); // points to first element of last cube  
   
    for (auto k= horizonLength - 1; k >= 0; --k, dst -= augR.n_rows * augA.n_rows * numModes) {
        dcube K_prev = K_k; // K_{k+1}
        
        dcube QAKA = augQ;
        dcube RBKB = augR;
        dcube BKA(augB.n_cols, augA.n_rows, numModes);
        
        for (uword i=0; i< numModes; ++i) {
            QAKA.slice(i) += symmatu(augA.slice(i).t() * K_prev.slice(i) * augA.slice(i)); // ensure symmetry
            RBKB.slice(i) += symmatu(augB.slice(i).t() * K_prev.slice(i) * augB.slice(i)); // ensure symmetry
            BKA.slice(i) = augB.slice(i).t() * K_prev.slice(i) * augA.slice(i);
        }
        
        dcube L(dst, augR.n_rows, augA.n_rows, numModes, false, true); //gains L_k for stage k
         
        for (uword j=0; j < numModes; ++j) {
            mat::const_row_iterator row = transitionMatrix.begin_row(j);
                        
            dmat P1 = zeros<dmat>(augA.n_rows, augA.n_cols);
            dmat P2 = zeros<dmat>(augB.n_cols, augA.n_rows);
            dmat P3 = zeros<dmat>(augR.n_rows, augR.n_cols);
            for (uword m = 0; m < numModes; ++m, ++row) {
                P1 += (*row) * QAKA.slice(m);
                P2 += (*row) * BKA.slice(m);
                P3 += (*row) * RBKB.slice(m);
            }
            K_k.slice(j) = P1 - symmatu(P2.t() * pinv(P3) * P2); // ensure symmetry
            L.slice(j) = -pinv(P3) * P2;            
        }      
    }
}

