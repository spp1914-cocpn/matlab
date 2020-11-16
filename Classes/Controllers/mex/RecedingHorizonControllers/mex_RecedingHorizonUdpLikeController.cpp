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

using namespace std;
using namespace arma;
using namespace RecedingHorizonUdpLikeController;

const double CONVERGENCE_DIFF = 1e-20;
const int MAX_NUM_ITERATIONS = 1000;

void computeGains(int stage, const Trajectory& trajectory, const Dynamics& dynamics, const Costate& costateEpsilon,  
       const dmat& JRJ, dmat& K, dmat& L, dmat& M) {
    
    const dcube& Xu = trajectory.getXu(stage);
    const dcube& Xl = trajectory.getXl(stage);
    const dmat& Se = dynamics.getSeAt(stage);
    const double* modeProbsCaPtr = dynamics.getModeProbsAt(stage);    
    
    // PlEpsilon and PuEpsilon are zero-ed first
    const dcube PuEpsilon = costateEpsilon.getPu();
    const dcube PlEpsilon = costateEpsilon.getPl();    

    const dmat SC = Se * dynamics.getAugC();
    
    dmat Psi = zeros<dmat>(Se.n_rows * PlEpsilon.n_rows, Se.n_cols * PlEpsilon.n_cols);
    dmat Ypsilon = zeros<dmat>(Psi.n_rows, Xl.n_cols * dynamics.getAugB(0).n_cols);
    dmat Pi = zeros<dmat>(Psi.n_rows, Xl.n_cols * PlEpsilon.n_cols);
    dmat Phi = symmatu(kron(Xl.slice(0), JRJ)); // JRJ is only non-zero for first mode
    dmat Sigma = zeros<dmat>(Phi.n_rows, Xl.n_cols * PlEpsilon.n_rows);
    dmat Lambda = zeros<dmat>(Sigma.n_cols, Xl.n_cols * PuEpsilon.n_cols);    
    
    dcolvec rho = zeros<dcolvec>(PlEpsilon.n_rows * SC.n_rows);
    dcolvec phi = zeros<dcolvec>(PlEpsilon.n_rows * Xl.n_cols);
    dcolvec gamma = zeros<dcolvec>(dynamics.getAugB(0).n_cols * Xl.n_cols);
    
    #ifdef _OPENMP
    // to use the reduction for Phi, a different initializer is required as it is initially not zero
    #pragma omp parallel for reduction(+:Psi,Ypsilon,Pi,Sigma,Phi,Lambda,rho,phi,gamma)
    #endif
    for (auto i=0; i < dynamics.getNumModes(); ++i) {
        const dmat noisePart = modeProbsCaPtr[i] * dynamics.getAugV();
               
        const dmat Xsum = Xl.slice(i) + Xu.slice(i);
        const dmat PlEpsilonB = PlEpsilon.slice(i) * dynamics.getAugB(i);
        const dmat PlEpsilonA = PlEpsilon.slice(i) * dynamics.getAugA(i);
        const dmat BPsum = dynamics.getAugB(i).t() * (PlEpsilon.slice(i) + PuEpsilon.slice(i));
        const dmat SCXl = SC * Xl.slice(i);
        
        Psi += symmatu(kron(Se * (noisePart + dynamics.getAugC() * Xsum * dynamics.getAugC().t()) * Se.t(), PlEpsilon.slice(i)));        
        Lambda += symmatu(kron(Xl.slice(i), PlEpsilon.slice(i)));
        Phi += symmatu(kron(Xl.slice(i), BPsum * dynamics.getAugB(i)));
        Ypsilon += kron(SCXl, PlEpsilonB);
        Sigma += kron(Xl.slice(i), PlEpsilonB.t());
        Pi += kron(SCXl, PlEpsilon.slice(i));
        // now come the vectors
        rho += vectorise(PlEpsilonA * Xsum * SC.t());
        phi += vectorise(PlEpsilonA * Xl.slice(i));
        gamma += vectorise(BPsum * dynamics.getAugA(i) * Xl.slice(i));        
    }
        
    const dmat A1 = join_horiz(join_horiz(Psi, -Ypsilon), Pi);
    const dmat A2 = join_horiz(join_horiz(-Ypsilon.t(), Phi), -Sigma);
    const dmat A3 = join_horiz(join_horiz(Pi.t(), -Sigma.t()), Lambda);
    const dmat A = symmatu(join_vert(join_vert(A1, A2), A3));
    const dcolvec b = join_vert(join_vert(-rho, gamma), -phi);
//     dcolvec gains;
//     if(solve(gains, A, -b)) {
//         K = reshape(gains.head(K.n_elem), K.n_rows, K.n_cols);
//         M = reshape(gains.tail(M.n_elem), M.n_rows, M.n_cols);
//         int startIdx = K.n_elem;
//         int endIdx = startIdx + L.n_elem-1;
//         L = reshape(gains.subvec(startIdx, endIdx), L.n_rows, L.n_cols);     
//     } else {
//         L.zeros();
//         K.zeros();
//         M.zeros();
//     }   

    dmat A_pinv(size(A));    
    if (pinv(A_pinv, A)) {
        // pseudo inverse could be obtained successfully
        const dcolvec gains = -A_pinv * b;
        K = reshape(gains.head(K.n_elem), K.n_rows, K.n_cols);
        M = reshape(gains.tail(M.n_elem), M.n_rows, M.n_cols);
        int startIdx = K.n_elem;
        int endIdx = startIdx + L.n_elem-1;
            
        L = reshape(gains.subvec(startIdx, endIdx), L.n_rows, L.n_cols);
    } else {
        // pseudo inverse could not be computed
        // fall back to zero gains
        K.zeros();
        M.zeros();
        L.zeros();
    }    
}

void mexFunction(int numOutputs, mxArray* outputArrays[], int numInputs, const mxArray* inputArrays[]) {
    
    // required params: augA, augB, augQ, augR, transitionMatrixCa, augW, augV, augC (sparse!), Se, JRJ, K, L, M, 
    // augStateCov, augStateSecondMoment, terminalAugQ, modeProbsCa
    auto augA = armaGetCubePr(inputArrays[0]);
    auto augB = armaGetCubePr(inputArrays[1]);
    auto augC = armaGetSparseMatrix(inputArrays[6]); // sparse
    auto augW = armaGetPr(inputArrays[4]);
    auto augV = armaGetPr(inputArrays[5]);
    auto transitionMatrixCa = armaGetPr(inputArrays[3]);
    auto modeProbsCa = armaGetPr(inputArrays[15]);
    auto Se = armaGetCubePr(inputArrays[7]);
    
    // augA, augB, augC, augW, augV, transitionMatrixCa, modeProbsCa, Se
    Dynamics dynamics(augA, augB, augC, augW, augV, transitionMatrixCa, modeProbsCa, Se);    
    
    const dcube augQ = armaGetCubePr(inputArrays[2]);
    const dmat JRJ = armaGetPr(inputArrays[8]);
    const dcube K_old = armaGetCubePr(inputArrays[9]);
    const dcube L_old = armaGetCubePr(inputArrays[10]);
    const dcube M_old = armaGetCubePr(inputArrays[11]);
    const dcube augStateCov_old = armaGetCubePr(inputArrays[12]);
    const dcube augStateSecondMoment_old = armaGetCubePr(inputArrays[13]);
    const dmat terminalAugQ = armaGetPr(inputArrays[14]);
    const int maxNumIterations = (numInputs > 16) ? armaGetScalar<int>(inputArrays[16]) : MAX_NUM_ITERATIONS;
    
    const int horizonLength = dynamics.getHorizonLength();
        
    size_t dims[3] = {static_cast<size_t>(K_old.n_rows), static_cast<size_t>(K_old.n_cols), static_cast<size_t>(horizonLength)};
    outputArrays[0] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // K
    dims[0] = static_cast<size_t>(L_old.n_rows);
    dims[1] = static_cast<size_t>(L_old.n_cols);
    outputArrays[1] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // L
    dims[0] = static_cast<size_t>(M_old.n_rows);
    dims[1] = static_cast<size_t>(M_old.n_cols);
    outputArrays[2] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // M
        
    dcube M_new(mxGetPr(outputArrays[2]), M_old.n_rows, M_old.n_cols, horizonLength, false, true);
    dcube L_new(mxGetPr(outputArrays[1]), L_old.n_rows, L_old.n_cols, horizonLength, false, true);
    dcube K_new(mxGetPr(outputArrays[0]), K_old.n_rows, K_old.n_cols, horizonLength, false, true);
    M_new = M_old;
    K_new = K_old;
    L_new = L_old;
    
    dcolvec oldCostToGo(horizonLength);
    oldCostToGo.fill(datum::inf);
    double oldCost = datum::inf;
    
    Trajectory predTrajectory(dynamics, augStateCov_old, augStateSecondMoment_old);
   
    int counter = 0;
    while (++counter < maxNumIterations) {
        // forward pass: obtain the second moments for the given sequence of controller gains
        predTrajectory.predictHorizon(K_new, L_new, M_new);
 
        Costate currentCostate = Costate::createTerminalCostate(dynamics, terminalAugQ);
         
        dmat tempK(K_old.n_rows, K_old.n_cols);
        dmat tempL(L_old.n_rows, L_old.n_cols);
        dmat tempM(M_old.n_rows, M_old.n_cols);
                
        // backward pass to compute the gains
        for (auto k= horizonLength - 1; k >= 0; --k) {
            // evaluate the epsilon operator            
            const Costate currentCostateEpsilon = currentCostate.evaluateEpsilonOperator();
                        
            computeGains(k, predTrajectory, dynamics, currentCostateEpsilon, JRJ, tempK, tempL, tempM);            
            
            const Costate tempCostate = currentCostateEpsilon.computeNew(tempK, tempL, tempM, augQ, JRJ);         

            double costToGo = predTrajectory.computeCostToGoForCostate(k, tempCostate);            
            if (costToGo < oldCostToGo.at(k)) {
                oldCostToGo(k) = costToGo;
                K_new.slice(k) = tempK;
                L_new.slice(k) = tempL;
                M_new.slice(k) = tempM;
                currentCostate = tempCostate;
            } else {
                // recompute the costate with the old gain                        
                currentCostate = currentCostateEpsilon.computeNew(K_new.slice(k), L_new.slice(k), M_new.slice(k),augQ, JRJ);
            }
        }
        if (abs(oldCost - oldCostToGo(0)) < CONVERGENCE_DIFF) {
            break;
        } else {
            oldCost = oldCostToGo(0);
        }        
    }    
    outputArrays[3] = mxCreateDoubleScalar(static_cast<double>(counter)); //numIterations
}

