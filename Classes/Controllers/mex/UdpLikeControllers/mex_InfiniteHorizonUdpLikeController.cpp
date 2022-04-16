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
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#pragma omp declare reduction(+ : dmat : omp_out += omp_in) initializer(omp_priv = zeros<dmat>(size(omp_orig)))
#endif

using namespace std;
using namespace arma;

const long MAX_ITER_NUM = 10000;
const long CONVERGENCE_DIFF = 1e-10;

void mexFunction(int numOutputs, mxArray* outputArrays[], int numInputs, const mxArray* inputArrays[]) {
        
    // required params: augA, augB, probC, expAugQ, expAugR, augW, augV, expAugA, expAugB, stationaryModeDist, augC
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const drowvec probC = armaGetPr(inputArrays[2]);    
    const dmat expAugQ = armaGetPr(inputArrays[3]);
    const dmat expAugR = armaGetPr(inputArrays[4]);    
    const dmat augW = armaGetPr(inputArrays[5]);
    const sp_dmat augV = armaGetSparseMatrix(inputArrays[6]); 
    const dmat expAugA = armaGetPr(inputArrays[7]);
    const dmat expAugB = armaGetPr(inputArrays[8]);
    const dcolvec stationaryModeDist = armaGetPr(inputArrays[9]);
    const dcube augC = armaGetCubePr(inputArrays[10]);        

    const uword dimInputSequence = expAugR.n_rows;
    const uword sequenceLength = augA.n_slices - 1;
    const uword dimAugState = expAugQ.n_rows;    
    const uword numCombinations = probC.n_elem;      
    const uword dimAugMeas = augC.n_rows;
        
    outputArrays[0] = mxCreateUninitNumericMatrix(static_cast<size_t>(dimInputSequence), static_cast<size_t>(dimAugState), mxDOUBLE_CLASS, mxREAL); // the gain L
    outputArrays[1] = mxCreateUninitNumericMatrix(static_cast<size_t>(dimAugState), static_cast<size_t>(dimAugMeas), mxDOUBLE_CLASS, mxREAL); // the gain K
    
    dmat L(mxGetPr(outputArrays[0]), dimInputSequence, dimAugState, false, true);
    dmat K(mxGetPr(outputArrays[1]), dimAugState, dimAugMeas, false, true);     
            
    dmat OverbarPsi = zeros<dmat>(dimAugState, dimAugState);
    dmat UnderbarPsi = zeros<dmat>(dimAugState, dimAugState);
    dmat OverbarLambda = zeros<dmat>(dimAugState, dimAugState);
    dmat UnderbarLambda = zeros<dmat>(dimAugState, dimAugState);    

    dmat OverbarPsi_old(size(OverbarPsi), fill::value(datum::inf));
    dmat UnderbarPsi_old(size(UnderbarPsi), fill::value(datum::inf));
         
    int counter = 0;
     
    while (counter < MAX_ITER_NUM
            && (sum(sum(abs(OverbarPsi_old - OverbarPsi))) > CONVERGENCE_DIFF || sum(sum(abs(UnderbarPsi_old - UnderbarPsi))) > CONVERGENCE_DIFF)) {
          
        // compute the next iterates
        const dmat K_old = K;
        const dmat L_old = L;
        
        OverbarPsi_old = OverbarPsi;
        UnderbarPsi_old = UnderbarPsi;
        const dmat OverbarLambda_old = OverbarLambda;
        const dmat UnderbarLambda_old = UnderbarLambda; 
        const dmat Lambda_sum = OverbarLambda_old + UnderbarLambda_old;
        
        const dmat BexpL = expAugB * L_old;
        const dmat part = expAugA - BexpL;
        const dmat part2 = symmatu(part * UnderbarPsi_old * part.t()); // should be symmetric        
        const dmat KVK = symmatu(K_old * augV * K_old.t()); // KVK'  should be symmetric
        const dmat LRL = symmatu(L_old.t() * expAugR * L_old); // L'RL should by symmetric
        
        UnderbarPsi = part2 + KVK;
        OverbarPsi = KVK + augW - part2;
        OverbarLambda = LRL + expAugQ - symmatu(part.t() * UnderbarLambda_old * part); // should be symmetric
        UnderbarLambda = LRL;
        
        // two separate terms required for L
        const dmat Bexp_Lambda = -expAugB.t() * UnderbarLambda_old;                              
        // compute the feedback gain L        
        dmat B_lambda_SumB = stationaryModeDist(0) * augB.slice(0).t() * Lambda_sum * augB.slice(0);
        dmat B_lambda_SumA = stationaryModeDist(0) * augB.slice(0).t() * Lambda_sum * augA.slice(0);
        for (uword i = 1; i <= sequenceLength; ++i) {
            B_lambda_SumB += stationaryModeDist(i) * augB.slice(i).t() * Lambda_sum * augB.slice(i);  
            B_lambda_SumA += stationaryModeDist(i) * augB.slice(i).t() * Lambda_sum * augA.slice(i);
        }
        const dmat L1 = Bexp_Lambda * expAugB + expAugR + symmatu(B_lambda_SumB);
        const dmat L2 = Bexp_Lambda * expAugA + B_lambda_SumA;        
        
        L = pinv(symmatu(L1)) * L2;            
                 
        // K_C for all i
        const dcube K_C = K_old * augC.each_slice();
        
        // two separate terms required for K
        dmat sumMat1 = zeros<dmat>(OverbarPsi_old.n_rows, augC.n_rows);      
        dmat sumMat2 = zeros<dmat>(augC.n_rows, augC.n_rows);        
        #ifdef _OPENMP    
        #pragma omp parallel for reduction(+:UnderbarLambda, sumMat1, sumMat2)
        #endif
        for (uword j = 0; j < numCombinations; ++j) {
            sumMat1 += probC(j) * OverbarPsi_old * augC.slice(j).t();
            sumMat2 += probC(j) * (augC.slice(j) * OverbarPsi_old * augC.slice(j).t());
            
            const dmat KC = K_C.slice(j);            
            const dmat partLambda = (expAugA - KC).t() * UnderbarLambda_old * (expAugA - KC) - (BexpL - KC).t() * UnderbarLambda_old * (BexpL - KC);            
            UnderbarLambda += probC(j) * symmatu(partLambda);            
        }        
        // expectation with respect to C        
        UnderbarPsi += symmatu(K_old * sumMat2 * K_old.t()); 
        
        //const dmat K1 = symmatu(augV + sumMat2);
        //const dmat K2 = expAugA * sumMat1;       
        K = (expAugA * sumMat1) * pinv(symmatu(augV + sumMat2));            
          
        #ifdef _OPENMP    
        #pragma omp parallel for reduction(+:UnderbarLambda, OverbarPsi, OverbarLambda)
        #endif
        for (uword i = 0; i <= sequenceLength; ++i) {           
            const dmat currAugA = augA.slice(i);
            const dmat currAugB = augB.slice(i);
            const dmat currBL = currAugB * L_old;            
                        
            dmat innerExpectationOverbarPsi = zeros<dmat>(currAugA.n_rows, currAugA.n_rows);
            dmat innerExpectationUnderbarLambda = zeros<dmat>(currBL.n_cols, currBL.n_cols);
            for (uword j = 0; j < numCombinations; ++j) {
                const dmat A_KC = currAugA - K_C.slice(j);
                const dmat BL_KC = currBL - K_C.slice(j);                
                
                innerExpectationOverbarPsi += probC(j) * symmatu(A_KC * OverbarPsi_old * A_KC.t());
                innerExpectationUnderbarLambda += probC(j) * symmatu(BL_KC.t() * UnderbarLambda_old * BL_KC);
            }           
            
            OverbarPsi += stationaryModeDist(i) * (innerExpectationOverbarPsi + symmatu((currAugA - currBL) * UnderbarPsi_old * (currAugA - currBL).t()));
            UnderbarLambda += stationaryModeDist(i) * (innerExpectationUnderbarLambda + symmatu(currBL.t() * OverbarLambda_old * currBL));
            OverbarLambda += stationaryModeDist(i) * symmatu((currAugA - currBL).t() * Lambda_sum * (currAugA - currBL));            
        }
        // end compute next iterates
        counter++;
    }
    L = -L;
}
