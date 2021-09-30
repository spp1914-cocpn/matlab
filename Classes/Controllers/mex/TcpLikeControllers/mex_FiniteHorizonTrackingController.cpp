/* >> This file is part of CoCPN-Sim
*
*    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
*
*    Copyright (C) 2017-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // required params: augA, augB, augQ, augR, transitionMatrix, terminalK, horizonLength, Qref, refWeightings
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const dcube augQ = armaGetCubePr(inputArrays[2]);
    const dcube augR = armaGetCubePr(inputArrays[3]);
    const dmat transitionMatrix = armaGetPr(inputArrays[4]);
    const dmat terminalK = armaGetPr(inputArrays[5]);
    const int horizonLength = armaGetScalar<int>(inputArrays[6]);
    const dmat Qref = armaGetPr(inputArrays[7]);
    const dmat refWeightings = armaGetPr(inputArrays[8]);
    
    const uword numModes = augA.n_slices;
    
    size_t dimsL[4] = {augR.n_rows, augA.n_rows, numModes, static_cast<size_t>(horizonLength)};
    outputArrays[0] = mxCreateUninitNumericArray(4, dimsL, mxDOUBLE_CLASS, mxREAL); // the gains L
    size_t dims[3] = {augR.n_rows, numModes, static_cast<size_t>(horizonLength)};
    outputArrays[1] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // the feedforward terms
    if (outputArrays[0] == NULL || outputArrays[1] == NULL) {
        // something went wrong
        mexErrMsgIdAndTxt("FiniteHorizonTrackingController:mex_FiniteHorizonTrackingController", 
                "** Allocation of output array failed **");
        return;
    }
    
    // we only need K_{k+1} to compute K_k and L_k
    dcube K_k(size(augA), fill::none); 
    K_k.each_slice() = terminalK; // K_N is the same for all modes
    // the same for sigma: only sigma_{k+1} required
    dmat sigma_k(refWeightings.n_rows, numModes, fill::none);
    sigma_k.each_col() = refWeightings.tail_cols(1);
    
    double* dstL = mxGetPr(outputArrays[0]);
    dstL += (dimsL[3] - 1) * dimsL[0] * dimsL[1] * dimsL[2]; // points to first element of last cube  of L
    
    double* dstFeedforward = mxGetPr(outputArrays[1]);
    dstFeedforward += (dims[2] - 1) * dims[0] * dims[1]; // points to first element of last slice of feedforward
            
    for (auto k = horizonLength - 1; k >= 0; --k, 
            dstL -= dimsL[0] * dimsL[1] * dimsL[2], 
            dstFeedforward -= dims[0] * dims[1]) {
        dcube K_prev = K_k; // K_{k+1}
        dmat sigma_prev = sigma_k; 
        
        dcube QAKA = augQ;
        dcube RBKB = augR;
        dcube BKA(augB.n_cols, augA.n_cols, numModes, fill::none);
        dmat Asigma(augA.n_cols, numModes, fill::none);
        dmat Bsigma(augB.n_cols, numModes, fill::none);
                
        for (uword i=0; i< numModes; ++i) {
            QAKA.slice(i) += symmatu(augA.slice(i).t() * K_prev.slice(i) * augA.slice(i)); // ensure symmetry
            RBKB.slice(i) += symmatu(augB.slice(i).t() * K_prev.slice(i) * augB.slice(i)); // ensure symmetry
            BKA.slice(i) = augB.slice(i).t() * K_prev.slice(i) * augA.slice(i);
            Asigma.col(i) = augA.slice(i).t() * sigma_prev.unsafe_col(i);
            Bsigma.col(i) = augB.slice(i).t() * sigma_prev.unsafe_col(i);
        }

        dcube L(dstL, dimsL[0], dimsL[1], dimsL[2], false, true); //gains L_k for stage k
        dmat feedforward(dstFeedforward, dims[0], dims[1], false, true); // feedforward terms for stage k
        
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (uword j=0; j < numModes; ++j) {
            mat::const_row_iterator row = transitionMatrix.begin_row(j);
            
            dmat P1 = zeros<dmat>(QAKA.n_rows, QAKA.n_cols);
            dmat P2 = zeros<dmat>(BKA.n_rows, BKA.n_cols);
            dmat P3 = zeros<dmat>(RBKB.n_rows, RBKB.n_cols);
            dcolvec s1 = zeros<dcolvec>(Asigma.n_rows);
            dcolvec s2 = zeros<dcolvec>(Bsigma.n_rows);
            
            for (uword m = 0; m < numModes; ++m, ++row) {
                P1 += (*row) * QAKA.slice(m);
                P2 += (*row) * BKA.slice(m);
                P3 += (*row) * RBKB.slice(m);
                s1 += (*row) * Asigma.unsafe_col(m);
                s2 += (*row) * Bsigma.unsafe_col(m);
            }
            K_k.slice(j) = P1 - symmatu(P2.t() * pinv(P3) * P2); // ensure symmetry
            sigma_k.col(j) = refWeightings.col(k) + s1 - P2.t() * pinv(P3) * s2;
            L.slice(j) = -pinv(P3) * P2;
            feedforward.col(j) = pinv(P3) * s2;            
        }      
    }
}

