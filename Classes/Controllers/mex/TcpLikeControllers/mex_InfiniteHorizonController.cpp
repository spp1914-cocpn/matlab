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

const double MIN_ITERATION_DIFF = 1e-6;
const long MAX_ITERATION_DIFF = 1e11;

int computeSteadyStateControlCovariance(const dcube& augA, const dcube& augB, const dcube& augQ, const dcube& augR, 
        const dmat& transitionMatrix, dcube& Pout) {
    
    const int numModes = augA.n_slices;
    
    int counter = 0;
    //drowvec convergenceDiffMin(numModes);
    dcube convergenceDiffMin(1,1, numModes);
    convergenceDiffMin.fill(MAX_ITERATION_DIFF);
        
    dcube currentP = zeros<dcube>(size(augA));
    dcube previousP = currentP;
    dcube Pmin = currentP;
    
    int status;
    
    while (true) {
        /* variable to monitor difference between two iterations of currP:
        we use sum sum of the difference as currP is known that PreviousP >
        currP in the positive definit sense. Therefore the quadratic form
        x' currP x with x = [1 1 .... 1 1]'. currP is sum(sum(P)).
        Therefore, this is an indicator for the change in positive-definite
        sense. */
        //const drowvec convergenceDiff = sum(sum(currentP - previousP));
        const dcube convergenceDiff = sum(sum(currentP - previousP, 1), 0);
        if (++counter > numModes + 1) {
            // ------- Terminate Conditions ------------------
            /*
             * Terminate Condition 1: stable / converged:
             * this happens if currP does not change anymore between two
             * iterations. Instead of zero, we use a lower bound
             */
            if (abs(convergenceDiff).max() < MIN_ITERATION_DIFF) {
                Pout = currentP;
                status = 1;
                break;
            } else if (abs(convergenceDiff).max() > MAX_ITERATION_DIFF) {                
                /*
                 * Terminate Condition 2: unstable / diverged:
                 * this happens if the closed-loop system is unstable (controller
                 * not able to stabilize system) or we had numerical problems that
                 * prevented to stop at the converged currP. The cases are hard to
                 * distinguish so we return the best solution found and the cautionflag set.
                 */
                Pout = Pmin;
                status = 0;
                break;
            } else if (convergenceDiff.min() <  0) {                
                /*
                 * Terminate Condition 3: Numerical Problems
                 * This should not happen but if it does we had numerical problems.
                 * The best solution found will be returned.
                 */
                Pout = Pmin;
                status = -1;
                break;
            }
            
            if (all(vectorise(convergenceDiff < convergenceDiffMin))) {
                Pmin = currentP;
                convergenceDiffMin = convergenceDiff;
            }
        }
        
        // computation within loop
        previousP = currentP;
        currentP.zeros();        
                    
        dcube QAPA = augQ;
        dcube RBPB = augR;
        dcube APB(augA.n_rows, augB.n_cols, numModes);
        
        for (auto i=0; i < numModes; ++i) {            
            QAPA.slice(i) += symmatu(augA.slice(i).t() * previousP.slice(i) * augA.slice(i));
            RBPB.slice(i) += symmatu(augB.slice(i).t() * previousP.slice(i) * augB.slice(i));
            APB.slice(i) = augA.slice(i).t() * previousP.slice(i) * augB.slice(i);
        }
        
        #ifdef _OPENMP
        #pragma omp parallel for shared(currentP, transitionMatrix, QAPA, RBPB, APB)
        #endif
        for (auto j=0; j < numModes; ++j) {
            mat::const_row_iterator row = transitionMatrix.begin_row(j);
            
            dmat P1 = zeros<dmat>(augQ.n_rows, augQ.n_cols);
            dmat P2 = zeros<dmat>(APB.n_rows, APB.n_cols);
            dmat P3 = zeros<dmat>(augR.n_rows, augR.n_cols);
            for (auto k=0; k < numModes; ++k, ++row) {
                P1 += (*row) * QAPA.slice(k);
                P2 += (*row) * APB.slice(k);
                P3 += (*row) * RBPB.slice(k);
            }           
            currentP.slice(j) = P1 - symmatu(P2 * pinv(P3) * P2.t());            
        }
    }    
    return status;
}

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // required params: augA, augB, augQ, augR, transitionMatrix, A
    const mwSize* const augADims = mxGetDimensions(inputArrays[0]);
    const mwSize* const augBDims = mxGetDimensions(inputArrays[1]);
    
    const int numModes = augADims[2];
    
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const dcube augQ = armaGetCubePr(inputArrays[2]);
    const dcube augR = armaGetCubePr(inputArrays[3]);
    const dmat transitionMatrix = armaGetPr(inputArrays[4]);
    
    size_t dims[3] = {static_cast<size_t>(augR.n_rows), static_cast<size_t>(augA.n_rows), static_cast<size_t>(numModes)};
    outputArrays[0] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // the gains L
        
    dcube L(mxGetPr(outputArrays[0]), augBDims[1], augADims[0], numModes, false, true);
        
    if (numInputs == 6) {
        const dmat A = armaGetPr(inputArrays[5]);
        const cx_vec eigenvalues = eig_gen(A);
        const double rho = max(abs(eigenvalues)); // spectral radius of A
        if (transitionMatrix(numModes - 1, numModes - 1) * rho * rho > 1) {
            // we return 
            int status = -2;
            outputArrays[1] = mxCreateDoubleScalar(static_cast<double>(status));
            // don't forget to copy the status variable to matlab
            //armaSetPr(outputArrays[1], status);
            return;
        }        
    }
 
    dcube Pout = zeros<dcube>(size(augA));
    int status = computeSteadyStateControlCovariance(augA, augB, augQ, augR, transitionMatrix, Pout);
    // now we can compute the actual gain matrices
    dcube RBPB = augR;
    dcube BPA(augB.n_cols, augA.n_rows, numModes);

    for (auto i=0; i < numModes; ++i) {
        RBPB.slice(i) += symmatu(augB.slice(i).t() * Pout.slice(i) * augB.slice(i));
        BPA.slice(i) = augB.slice(i).t() * Pout.slice(i) * augA.slice(i);
    }

    for (auto j=0; j < numModes; ++j) {
        mat::const_row_iterator row = transitionMatrix.begin_row(j);

        dmat L1 = zeros<dmat>(augR.n_rows, augR.n_rows);
        dmat L2 = zeros<dmat>(BPA.n_rows, BPA.n_cols);
        for (auto k=0; k < numModes; ++k, ++row) {
            L1 += (*row) * RBPB.slice(k);
            L2 += (*row) * BPA.slice(k);                
        }           
        L.slice(j) = -pinv(L1) * L2;
    }
    // don't forget to copy the status variable to matlab
    outputArrays[1] = mxCreateDoubleScalar(static_cast<double>(status));
}

