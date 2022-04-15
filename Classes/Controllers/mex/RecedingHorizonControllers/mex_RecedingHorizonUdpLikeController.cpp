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

#include "include/Costate.h"
#include "include/Trajectory.h"
#include "include/Dynamics.h"
#include <armadillo>
#include "armaMex.hpp"
#ifdef _OPENMP
#include <omp.h>
#pragma omp declare reduction(+: dmat : omp_out += omp_in) initializer(omp_priv = zeros<dmat>(size(omp_orig)))
#pragma omp declare reduction(+: dcolvec : omp_out += omp_in) initializer(omp_priv = omp_orig)
#endif

#ifdef __cplusplus
extern "C" bool utIsInterruptPending();
#else
extern bool utIsInterruptPending();
#endif

using namespace std;
using namespace arma;
using namespace RecedingHorizonUdpLikeController;

const double CONVERGENCE_DIFF = 1e-20;
const int MAX_NUM_ITERATIONS = 5;

const solve_opts::opts SOLVE_OPTS = solve_opts::refine + solve_opts::equilibrate + solve_opts::allow_ugly;

void computeGains(int stage, const Trajectory& trajectory, const Dynamics& dynamics, const Costate& costateEpsilon,  
       const dmat& JRJ, dmat& K, dmat& L) {
    
    const dcube& Xu = trajectory.getXu(stage);
    const dcube& Xl = trajectory.getXl(stage);
    const dmat& Se = dynamics.getSeAt(stage);
    const dmat& Aexp = dynamics.getAexpAt(stage);
    const dmat& Bexp = dynamics.getBexpAt(stage);
    const double* modeProbsCaPtr = dynamics.getModeProbsAt(stage);    
    
    // PlEpsilon and PuEpsilon are zero-ed first
    const dcube PuEpsilon = costateEpsilon.getPu();
    const dcube PlEpsilon = costateEpsilon.getPl();    

    const dmat CS = (Se * dynamics.getAugC()).t();
    
    // for K_k
    dmat Psi = zeros<dmat>(Se.n_rows * PlEpsilon.n_rows, Se.n_cols * PlEpsilon.n_cols);    
    // for L_k
    dmat Phi = symmatu(kron(Xl.slice(0), JRJ)); // JRJ is only non-zero for first mode
    
    // for K_k
    // negative
    dcolvec rho = zeros<dcolvec>(PlEpsilon.n_rows * CS.n_cols);    
    // for L_k
    dcolvec gamma = zeros<dcolvec>(dynamics.getAugB(0).n_cols * Xl.n_cols);
    
    #ifdef _OPENMP
    // to use the reduction for Phi, a different initializer is required as it is initially not zero
    #pragma omp parallel for reduction(+:Psi,Phi,rho,gamma)
    #endif
    for (auto i=0; i < dynamics.getNumModes(); ++i) {
        const dmat noisePart = modeProbsCaPtr[i] * dynamics.getAugV();  
        
        // quadratic term for K_k
        Psi += symmatu(kron(Se * (noisePart + dynamics.getAugC() * Xu.slice(i) * dynamics.getAugC().t()) * Se.t(), PlEpsilon.slice(i)));        
        // quadratic term for L_k
        const dmat diffBTranspose = (dynamics.getAugB(i) - Bexp).t();
        const dmat partPhi = symmatu(dynamics.getAugB(i).t() * PuEpsilon.slice(i) * dynamics.getAugB(i) 
            + diffBTranspose * PlEpsilon.slice(i) * diffBTranspose.t());
        Phi += kron(Xl.slice(i), partPhi);
        
        // now come the vectors
        // this is -rho
        rho += vectorise(PlEpsilon.slice(i) * dynamics.getAugA(i) * Xu.slice(i) * CS);
        const dmat partGamma = dynamics.getAugB(i).t() * PuEpsilon.slice(i) * dynamics.getAugA(i) 
            + diffBTranspose * PlEpsilon.slice(i) * (dynamics.getAugA(i) - Aexp);
        gamma += vectorise(partGamma * Xl.slice(i));        
    }

    // we can solve independently for K and L
    // first, for K
    dmat psiInv(size(Psi), fill::none);
    if (rho.is_zero()) {
        // shortcut: rho is zero, so return zero solution
        K.zeros();
    } else if (pinv(psiInv, Psi)) {
        // pseudo inverse could be obtained successfully
        const dcolvec gainK = psiInv * rho; // rho is negative
        K = reshape(gainK, K.n_rows, K.n_cols);
    } else {
        // try to use solve directly
        dcolvec gainK;
        if (solve(gainK, Psi, rho, SOLVE_OPTS) && gainK.is_finite()) {
            K = reshape(gainK, K.n_rows, K.n_cols);
        }
        else {
            mexWarnMsgIdAndTxt("RecedingHorizonController:mex_RecedingHorizonController" ,
                "** K could not be updated: both pinv and solve failed **");
        }
    }

    if (gamma.is_zero()) {
        // shortcut: gamma is zero, so return zero solution
        L.zeros();
        return;
    }    
    // now solve for L
    dmat phiInv(size(Phi), fill::none);
    if (pinv(phiInv, Phi)) {
        // pseudo inverse could be obtained successfully
        const dcolvec gainL = -phiInv * gamma;
        L = reshape(gainL, L.n_rows, L.n_cols);        
    } else {
        // try to use solve directly
        dcolvec gainL;
        if (solve(gainL, Phi, -gamma, SOLVE_OPTS)  && gainL.is_finite()) {
            L = reshape(gainL, L.n_rows, L.n_cols);  
        }
        else {
            mexWarnMsgIdAndTxt("RecedingHorizonController:mex_RecedingHorizonController" ,
                "** L could not be updated: both pinv and solve failed **");
        }
    }
}

void mexFunction(int numOutputs, mxArray* outputArrays[], int numInputs, const mxArray* inputArrays[]) {
    
    // required params: augA, augB, augQ, transitionMatrixCa, augW, augV, augC (sparse!), Se, JRJ, K, L, 
    // augStateCov, augStateSecondMoment, terminalAugQ, modeProbsCa
    dcube augA = armaGetCubePr(inputArrays[0]);
    dcube augB = armaGetCubePr(inputArrays[1]);
    sp_dmat augC = armaGetSparseMatrix(inputArrays[6]); // sparse
    dmat augW = armaGetPr(inputArrays[4]);
    dmat augV = armaGetPr(inputArrays[5]);
    dmat transitionMatrixCa = armaGetPr(inputArrays[3]);
    dcolvec modeProbsCa = armaGetPr(inputArrays[14]);
    dcube Se = armaGetCubePr(inputArrays[7]);
   
    // augA, augB, augC, augW, augV, transitionMatrixCa, modeProbsCa, Se
    Dynamics dynamics(augA, augB, augC, augW, augV, transitionMatrixCa, modeProbsCa, Se);    
    
    const dcube augQ = armaGetCubePr(inputArrays[2]);
    const dmat JRJ = armaGetPr(inputArrays[8]);
    const dcube K_old = armaGetCubePr(inputArrays[9]);
    const dcube L_old = armaGetCubePr(inputArrays[10]);
    const dcube augStateCov_old = armaGetCubePr(inputArrays[11]);
    const dcube augStateSecondMoment_old = armaGetCubePr(inputArrays[12]);
    const dmat terminalAugQ = armaGetPr(inputArrays[13]);
    const int maxNumIterations = (numInputs > 15) ? armaGetScalar<int>(inputArrays[15]) : MAX_NUM_ITERATIONS;
    
    const int horizonLength = dynamics.getHorizonLength();
        
    size_t dims[3] = {static_cast<size_t>(K_old.n_rows), static_cast<size_t>(K_old.n_cols), static_cast<size_t>(horizonLength)};
    outputArrays[0] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // K
    dims[0] = static_cast<size_t>(L_old.n_rows);
    dims[1] = static_cast<size_t>(L_old.n_cols);
    outputArrays[1] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // L    
    
    dcube L_new(mxGetPr(outputArrays[1]), L_old.n_rows, L_old.n_cols, horizonLength, false, true);
    dcube K_new(mxGetPr(outputArrays[0]), K_old.n_rows, K_old.n_cols, horizonLength, false, true);
    K_new = K_old;
    L_new = L_old;
    
    dcolvec costToGo(horizonLength);
    costToGo.fill(datum::inf);
    double cost = datum::inf;
    
    Trajectory predTrajectory(dynamics, augStateCov_old, augStateSecondMoment_old);
   
    int counter = 0;
    while (counter < maxNumIterations) {
        
        counter++;
        // forward pass: obtain the second moments for the given sequence of controller gains
        predTrajectory.predictHorizon(K_new, L_new);
 
        Costate currentCostate = Costate::createTerminalCostate(dynamics, terminalAugQ);
         
        dmat tempK(K_old.n_rows, K_old.n_cols);
        dmat tempL(L_old.n_rows, L_old.n_cols);
        
        // backward pass to compute the gains
        for (auto k= horizonLength - 1; k >= 0; --k) {
            // evaluate the epsilon operator            
            const Costate currentCostateEpsilon = currentCostate.evaluateEpsilonOperator();
                       
            computeGains(k, predTrajectory, dynamics, currentCostateEpsilon, JRJ, tempK, tempL);            
            
            const Costate tempCostate = currentCostateEpsilon.computeNew(tempK, tempL, augQ, JRJ);         

            double currCostToGo = predTrajectory.computeCostToGoForCostate(k, tempCostate);
            if (currCostToGo < costToGo(k)) {
                costToGo(k) = currCostToGo;
                K_new.slice(k) = tempK;
                L_new.slice(k) = tempL;      
                currentCostate = tempCostate;
            } else {
                // recompute the costate with the old gain                        
                currentCostate = currentCostateEpsilon.computeNew(K_new.slice(k), L_new.slice(k), augQ, JRJ);
            }            
        }
        
        const double terminate = abs(cost - costToGo(0)) < CONVERGENCE_DIFF;
        cost = costToGo(0);
        if (terminate) {
            break;
        }
         
        if (utIsInterruptPending()) {        /* check for a Ctrl-C event */
            mexPrintf("Ctrl-C Detected. END\n\n");
            break;
        }
    }
    mexPrintf("** RecedingHorizonUdpLikeController: %d iterations, current costs (bound): %f **\n", counter, cost);
    
    outputArrays[2] = mxCreateDoubleScalar(static_cast<double>(counter)); //numIterations
    outputArrays[3] = mxCreateUninitNumericMatrix(static_cast<size_t>(augA.n_rows),
            static_cast<size_t>(augA.n_cols), mxDOUBLE_CLASS, mxREAL); // expAugA[0]
    armaSetPr(outputArrays[3], dynamics.getAexpAt(0));
    outputArrays[4] = mxCreateUninitNumericMatrix(static_cast<size_t>(augB.n_rows),
            static_cast<size_t>(augB.n_cols), mxDOUBLE_CLASS, mxREAL); // expAugB[0]
    armaSetPr(outputArrays[4], dynamics.getBexpAt(0));
}

