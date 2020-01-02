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
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // required params: augA, augB, augQ, augR, transitionMatrix, terminalK, terminalP, 
    // inputWeightings, stateWeightings, horizonLength
            
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const dcube augQ = armaGetCubePr(inputArrays[2]);
    const dcube augR = armaGetCubePr(inputArrays[3]);
    const dmat transitionMatrix = armaGetPr(inputArrays[4]);
    const dmat terminalK = armaGetPr(inputArrays[5]);
    const dmat terminalP = armaGetPr(inputArrays[6]);
    
    dcube inputWeightings, stateWeightings;
    const int ndims = mxGetNumberOfDimensions(inputArrays[7]);
    const mwSize* inputDims = mxGetDimensions(inputArrays[7]);
    const mwSize* stateDims = mxGetDimensions(inputArrays[8]);
    if (ndims == 3){
        inputWeightings = dcube(mxGetPr(inputArrays[7]), inputDims[0], inputDims[1], inputDims[2], false, true);
        stateWeightings = dcube(mxGetPr(inputArrays[8]), stateDims[0], stateDims[1], stateDims[2], false, true);
    } else {
        // only one constraint present
        inputWeightings = dcube(mxGetPr(inputArrays[7]), inputDims[0], inputDims[1], 1, false, true);
        stateWeightings = dcube(mxGetPr(inputArrays[8]), stateDims[0], stateDims[1], 1, false, true);
    }   
        
    const int horizonLength = armaGetScalar<int>(inputArrays[9]);
    
    const uword numModes = augA.n_slices;
    const int sequenceLength = numModes - 1;
    const uword numConstraints = stateWeightings.n_slices;
    const uword dimPlantInput = inputWeightings.n_rows;
    const uword dimPlantState = stateWeightings.n_rows;
    const uword dimAugState = augA.n_rows;
            
    size_t dims[4] = {augR.n_rows, dimAugState, numModes, static_cast<size_t>(horizonLength)};
    outputArrays[0] = mxCreateUninitNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL); // the gains L
    // this one requires a lot of memory
    dims[1] = numConstraints;
    outputArrays[1] = mxCreateUninitNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL); // the feedforward terms
    size_t dimsP[3] = {terminalP.n_rows, terminalP.n_cols, numModes};
    outputArrays[2] = mxCreateUninitNumericArray(3, dimsP, mxDOUBLE_CLASS, mxREAL); // P_0
    size_t dimsS[3] = {numConstraints, numConstraints, numModes};
    outputArrays[3] = mxCreateNumericArray(3, dimsS, mxDOUBLE_CLASS, mxREAL); // S_0, initialied to zero
    
    if (outputArrays[0] == NULL || outputArrays[1] == NULL  || outputArrays[2] == NULL || outputArrays[3] == NULL) {
        // something went wrong
        mexErrMsgIdAndTxt("FiniteHorizonController:mex_FiniteHorizonControllerWithConstraints", 
                "** Allocation of output arrays failed **");
        return;
    }        
    
    // we only need S_0 in the end for the computation of the multiplier, so just store S_k and S_{k+1} per mode
    dcube S_k(mxGetPr(outputArrays[3]), numConstraints, numConstraints, numModes, false, true);
    // we only need P_0 in the end for the computation of the multiplier, so just store P_k and P_{k+1} per mode
    dcube P_k(mxGetPr(outputArrays[2]), terminalP.n_rows, terminalP.n_cols, numModes, false, true);
    P_k.each_slice() = terminalP; // P_tilde_N in the paper is the same for all modes
    
    // we only need K_{k+1} to compute K_k and L_k
    dcube K_k(size(augA)); 
    K_k.each_slice() = terminalK; // K_N is the same for all modes
        
    double* dstL = mxGetPr(outputArrays[0]);
    dstL += ((horizonLength-1) * augR.n_rows * dimAugState * numModes); // points to first element of last cube of L
   
    double* dstFeedforward = mxGetPr(outputArrays[1]);
    dstFeedforward += ((horizonLength-1) * augR.n_rows * numConstraints * numModes); // points to first element of last cube  
    
    const uvec sums = cumsum(regspace<uvec>(1, 1, sequenceLength - 1)) * dimPlantInput; // explicitly pass delta = 1 to be compatible with existing Matlab implementation (empty vector, if sequenceLength -1 > 1)
    uvec startIdx(sums.n_elem);
    int i = 2;
    startIdx.imbue([&] () {return dimAugState-sums(sequenceLength-(i++));});
    const uvec endIdx = startIdx + dimPlantInput - 1;
    
    for (auto k= horizonLength - 1; k >= 0; --k, dstL -= augR.n_rows * dimAugState * numModes, dstFeedforward -=augR.n_rows * numConstraints * numModes) {
        const dcube K_prev = K_k; // K_{k+1}
        const dcube S_prev = S_k;
        const dcube P_prev = P_k;
        
        const dmat inputWeights = inputWeightings.col(k);
        const dmat stateWeights = stateWeightings.col(k);        
        
        dcube Q_tilde = zeros<dcube>(augA.n_rows, numConstraints, numModes);                
        Q_tilde.each_slice([&stateWeights, &dimPlantState](dmat& s) {s.head_rows(dimPlantState) = stateWeights;}, true);
        
        // Q_tilde (per mode) depends on H matrix which is
        // zero for first and last mode
        // inline code (with index shift) from Utility.createAugmentedLinearIntegralConstraints
        // to compute the mode-dependent Q_tilde
        for (uword i=1; i < numModes -1; ++i) {           
            Q_tilde(span(startIdx(i-1), endIdx(i-1)), span::all, span(i)) = inputWeights;            
        }
        
        dcube QAKA = augQ;
        dcube RBKB = augR;
        dcube BKA(augB.n_cols, dimAugState, numModes);
        dcube APQ_tilde = Q_tilde;
        dcube BPR_tilde(augB.n_cols, numConstraints, numModes); 
        
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (uword i=0; i< numModes; ++i) {
            QAKA.slice(i) += symmatu(augA.slice(i).t() * K_prev.slice(i) * augA.slice(i)); // ensure symmetry
            RBKB.slice(i) += symmatu(augB.slice(i).t() * K_prev.slice(i) * augB.slice(i)); // ensure symmetry
            BKA.slice(i) = augB.slice(i).t() * K_prev.slice(i) * augA.slice(i);
            APQ_tilde.slice(i) += augA.slice(i).t() * P_prev.slice(i);
            BPR_tilde.slice(i) = augB.slice(i).t() * P_prev.slice(i);
        }
        //R_tilde (augmented input constraint weightings) only nonzero for the first mode
        BPR_tilde.slice(0).head_rows(dimPlantInput) += inputWeights;
        
        dcube L(dstL, augR.n_rows, dimAugState, numModes, false, true); //gains L_k for stage k
        dcube feedforward(dstFeedforward, augR.n_rows, numConstraints, numModes, false, true); // feedforward terms for stage k
        
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (uword j=0; j < numModes; ++j) {
            mat::const_row_iterator row = transitionMatrix.begin_row(j);
                        
            dmat P1 = zeros<dmat>(dimAugState, dimAugState);
            dmat P2 = zeros<dmat>(augB.n_cols, dimAugState);
            dmat P3 = zeros<dmat>(augR.n_rows, augR.n_cols);
            dmat P4 = zeros<dmat>(APQ_tilde.n_rows, APQ_tilde.n_cols);
            dmat P5 = zeros<dmat>(BPR_tilde.n_rows, BPR_tilde.n_cols);
            dmat S1 = zeros<dmat>(numConstraints, numConstraints);
            for (uword m = 0; m < numModes; ++m, ++row) {
                P1 += (*row) * QAKA.slice(m);
                P2 += (*row) * BKA.slice(m);
                P3 += (*row) * RBKB.slice(m);
                P4 += (*row) * APQ_tilde.slice(m);
                P5 += (*row) * BPR_tilde.slice(m);
                
                S1 += (*row) * S_prev.slice(m);
            }
            K_k.slice(j) = P1 - symmatu(P2.t() * pinv(P3) * P2); // ensure symmetry
            P_k.slice(j) = P4 - P2.t() * pinv(P3) * P5;
            L.slice(j) = -pinv(P3) * P2;
            feedforward.slice(j) = -pinv(P3) * P5;
            S_k.slice(j) = S1 - symmatu(P5.t() * pinv(P3) * P5);
        }      
    }
}

