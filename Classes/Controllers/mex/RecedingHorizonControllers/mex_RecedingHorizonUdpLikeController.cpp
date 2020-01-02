/* >> This file is part of CoCPN-Sim
*
*    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
*
*    Copyright (C) 2018-2019  Florian Rosenthal <florian.rosenthal@kit.edu>
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
#pragma omp declare reduction(+: dmat : omp_out += omp_in) initializer(omp_priv = zeros<dmat>(size(omp_orig)))
#pragma omp declare reduction(+: dcolvec : omp_out += omp_in) initializer(omp_priv = omp_orig)
#endif

using namespace std;
using namespace arma;

struct costate {
    dcube Pu;
    dcube Pl;
    dcolvec omega;
};

const double CONVERGENCE_DIFF = 1e-10;
const int MAX_NUM_ITERATIONS = 1000;


void computeGains(const dcube& Xu, const dcube& Xl, const costate& costateEpsilon, const double* modeProbsCaPtr, 
        const dmat& transitionMatrixCa, const dcube& augA, const dcube& augB, const sp_dmat& augC, 
        const dmat& augV, const dmat& Se, const dmat& JRJ, 
        dmat& K, dmat& L, dmat& M) {
        
    // PlEpsilon and PuEpsilon are zero-ed first
    const dcube PuEpsilon = costateEpsilon.Pu;
    const dcube PlEpsilon = costateEpsilon.Pl;
    
    const int numCaModes = transitionMatrixCa.n_rows;
    const dmat SC = Se * augC;
    
    dmat Psi = zeros<dmat>(Se.n_rows * PlEpsilon.n_rows, Se.n_cols * PlEpsilon.n_cols);
    dmat Ypsilon = zeros<dmat>(Psi.n_rows, Xl.n_cols * augB.n_cols);
    dmat Pi = zeros<dmat>(Psi.n_rows, Xl.n_cols * PlEpsilon.n_cols);
    dmat Phi = symmatu(kron(Xl.slice(0), JRJ)); // JRJ is only non-zero for first mode
    dmat Sigma = zeros<dmat>(Phi.n_rows, Xl.n_cols * PlEpsilon.n_rows);
    dmat Lambda = zeros<dmat>(Sigma.n_cols, Xl.n_cols * PuEpsilon.n_cols);    
    
    dcolvec rho = zeros<dcolvec>(PlEpsilon.n_rows * SC.n_rows);
    dcolvec phi = zeros<dcolvec>(PlEpsilon.n_rows * Xl.n_cols);
    dcolvec gamma = zeros<dcolvec>(augB.n_cols * Xl.n_cols);
    
    #ifdef _OPENMP
    // to use the reduction for Phi, a different initializer is required as it is initially not zero
    #pragma omp parallel for reduction(+:Psi,Ypsilon,Pi,Sigma,Phi,Lambda,rho,phi,gamma) \
        shared(transitionMatrixCa, modeProbsCaPtr, Xu, Xl)
    #endif
    for (auto i=0; i < numCaModes; ++i) {
        const dmat noisePart = modeProbsCaPtr[i] * augV;
               
        const dmat Xsum = Xl.slice(i) + Xu.slice(i);
        const dmat PlEpsilonB = PlEpsilon.slice(i) * augB.slice(i);
        const dmat PlEpsilonA = PlEpsilon.slice(i) * augA.slice(i);
        const dmat BPsum = augB.slice(i).t() * (PlEpsilon.slice(i) + PuEpsilon.slice(i));
        const dmat SCXl = SC * Xl.slice(i);
        
        Psi += symmatu(kron(Se * (noisePart + augC * Xsum * augC.t()) * Se.t(), PlEpsilon.slice(i)));        
        Lambda += symmatu(kron(Xl.slice(i), PlEpsilon.slice(i)));
        Phi += symmatu(kron(Xl.slice(i), BPsum * augB.slice(i)));
        Ypsilon += kron(SCXl, PlEpsilonB);
        Sigma += kron(Xl.slice(i), PlEpsilonB.t());
        Pi += kron(SCXl, PlEpsilon.slice(i));
        // now come the vectors
        rho += vectorise(PlEpsilonA * Xsum * SC.t());
        phi += vectorise(PlEpsilonA * Xl.slice(i));
        gamma += vectorise(BPsum * augA.slice(i) * Xl.slice(i));        
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
        // fal back to zero gains
        K.zeros();
        M.zeros();
        L.zeros();
    }    
}

const costate computeCostate(const costate& costateEpsilon,
        const dcube& augA, const dcube& augB, const dmat& K, const dmat& L, const dmat& M, 
        const dmat& augV, const dmat& augW, const dmat& Se, const sp_dmat& augC, const dmat& transitionMatrixCa,
        const dcube& augQ, const dmat& JRJ) {
    
    const int numCaModes = transitionMatrixCa.n_rows;
        
    costate result = {augQ, zeros<dcube>(size(augQ)), transitionMatrixCa * costateEpsilon.omega};
    result.Pu.slice(0) += symmatu(L.t() * JRJ * L);
    result.Pl.slice(0) += symmatu(L.t() * JRJ * L);
            
    const dmat KS = K * Se;
    const dmat KSC = KS * augC;
    const dmat E_tilde = symmatu(KS * augV * KS.t()); // E_tilde in the paper            
    
    #ifdef _OPENMP
    #pragma omp parallel for shared(result, costateEpsilon)
    #endif
    for (int i = 0; i < numCaModes; ++i) {
        const dmat currentPlEpsilon = costateEpsilon.Pl.slice(i);
        const dmat currentPuEpsilon = costateEpsilon.Pu.slice(i);
        const dmat BL = augB.slice(i) * L;
        const dmat U_tilde = augA.slice(i) + BL; // U_tilde in the paper
        const dmat O_tilde = M - BL; // O_tilde in the paper
        const dmat D_tilde = U_tilde - M - KSC; //D_tilde in the paper
        
        result.omega(i) += trace(currentPlEpsilon * E_tilde + (currentPlEpsilon + currentPuEpsilon) * augW);
        result.Pu.slice(i) += symmatu(U_tilde.t() * currentPuEpsilon * U_tilde + D_tilde.t() * currentPlEpsilon * D_tilde);
        result.Pl.slice(i) += symmatu(BL.t() * currentPuEpsilon * BL + O_tilde.t() * currentPlEpsilon * O_tilde);        
    }
    return result;
}

void predictXHorizon(const dcube& augA, const dcube& augB, const dcube& K, const dcube& L, const dcube& M, 
        const dmat& augV, const dmat& augW, const dcube& Se, const sp_dmat& augC, const dmat& modeProbsCa, 
        const dmat& transitionMatrixCa, field<cube>& Xu, field<cube>& Xl) {
    
    const int numCaModes = modeProbsCa.n_rows;
    const int horizonLength = Se.n_slices;
    //first element (cube) of Xu, Xl is already filled
        
    for (auto k= 0; k < horizonLength; ++k) {
        const dmat KS = K.slice(k) * Se.slice(k);
        const dmat KSC = KS * augC;
        const dmat E_tilde = symmatu(KS * augV * KS.t()); // E_tilde in the paper
        const dmat MKSC = M.slice(k) + KSC;
        
        const dcube currXu = Xu.at(k);
        const dcube currXl = Xl.at(k);
        const double* modeCol = modeProbsCa.colptr(k);
        dcube firstPart(E_tilde.n_rows, E_tilde.n_cols, numCaModes);
        dcube secondPart(E_tilde.n_rows, E_tilde.n_cols, numCaModes);
        
        #ifdef _OPENMP
        #pragma omp parallel for shared(firstPart, secondPart, augW)
        #endif
        for (auto i = 0; i < numCaModes; ++i) {
            const dmat AKSC = augA.slice(i) - KSC;
            const dmat D_tilde = AKSC - M.slice(k) + augB.slice(i) * L.slice(k); // D_tilde in the paper
            firstPart.slice(i) = symmatu(AKSC * currXu.slice(i) * AKSC.t() + D_tilde * currXl.slice(i) * D_tilde.t())
                            + modeCol[i] * (E_tilde + augW);
            secondPart.slice(i) = symmatu(MKSC * currXl.slice(i) * MKSC.t() + KSC * currXu.slice(i) * KSC.t())
                            + modeCol[i] * E_tilde;
        }
        
        Xu.at(k+1).zeros();
        Xl.at(k+1).zeros();
        
        #ifdef _OPENMP        
        #pragma omp parallel for
        #endif
        for (auto j=0; j < numCaModes; ++j) { // new mode            
            const double* col = transitionMatrixCa.colptr(j);
            for (auto i=0; i < numCaModes; ++i) { // old mode
                Xu.at(k+1).slice(j) += col[i] * firstPart.slice(i);
                Xl.at(k+1).slice(j) += col[i] * secondPart.slice(i);
            }
        }
    }    
}

const costate evaluateEpsilonOperator(const costate& forCostate, const dmat& transitionMatrixCa) {
    
    const int numCaModes = transitionMatrixCa.n_rows;
        
    costate result = {zeros<dcube>(size(forCostate.Pu)), zeros<dcube>(size(forCostate.Pl)), 
        transitionMatrixCa * forCostate.omega};
        
    for (auto i=0; i < numCaModes; ++i) {
        mat::const_row_iterator row = transitionMatrixCa.begin_row(i);
        for (auto j=0; j < numCaModes; ++j, ++row) {
            result.Pu.slice(i) += (*row) * forCostate.Pu.slice(j);
            result.Pl.slice(i) += (*row) * forCostate.Pl.slice(j);            
        }  
    }
    return result;
}

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // required params: augA, augB, augQ, augR, transitionMatrixCa, augW, augV, augC, Se, JRJ, K, L, M, 
    // augStateCov, augStateSecondMoment, terminalAugQ, modeProbsCa
        
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const dcube augQ = armaGetCubePr(inputArrays[2]);
    const dmat transitionMatrixCa = armaGetPr(inputArrays[3]);
    const dmat augW = armaGetPr(inputArrays[4]);
    const dmat augV = armaGetPr(inputArrays[5]);
    const sp_dmat augC = armaGetSparseMatrix(inputArrays[6]); // this one is sparse
    const dcube Se = armaGetCubePr(inputArrays[7]);
    const dmat JRJ = armaGetPr(inputArrays[8]);
    const dcube K_old = armaGetCubePr(inputArrays[9]);
    const dcube L_old = armaGetCubePr(inputArrays[10]);
    const dcube M_old = armaGetCubePr(inputArrays[11]);
    const dcube augStateCov_old = armaGetCubePr(inputArrays[12]);
    const dcube augStateSecondMoment_old = armaGetCubePr(inputArrays[13]);
    const dmat terminalAugQ = armaGetPr(inputArrays[14]);
    const dmat modeProbsCa = armaGetPr(inputArrays[15]);
    
    const int horizonLength = Se.n_slices;
    const int numCaModes = modeProbsCa.n_rows;
    
    size_t dims[3] = {static_cast<size_t>(K_old.n_rows), static_cast<size_t>(K_old.n_cols), static_cast<size_t>(horizonLength)};
    outputArrays[0] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // K
    dims[0]  = static_cast<size_t>(L_old.n_rows);
    dims[1] =  static_cast<size_t>(L_old.n_cols);
    outputArrays[1] = mxCreateUninitNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL); // L
    dims[0]  = static_cast<size_t>(M_old.n_rows);
    dims[1] =  static_cast<size_t>(M_old.n_cols);
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
    
    field<dcube> Xu(horizonLength + 1);
    field<dcube> Xl(horizonLength + 1);
    Xu.for_each([&augStateCov_old](dcube& c) {c.copy_size(augStateCov_old);});
    Xl.for_each([&augStateSecondMoment_old](dcube& c) {c.copy_size(augStateSecondMoment_old);});
    Xu.at(0) = augStateCov_old;
    Xl.at(0) = augStateSecondMoment_old;
           
    int counter = 0;
    while (++counter < MAX_NUM_ITERATIONS) {
        // forward pass: obtain the second moments for the given sequence of controller gains
        predictXHorizon(augA, augB, K_new, L_new, M_new, augV, augW, Se, augC, 
                modeProbsCa, transitionMatrixCa, Xu, Xl);
        
        costate currentCostate = {dcube(size(augA)), zeros<dcube>(size(augA)), zeros<dcolvec>(numCaModes)};
        currentCostate.Pu.each_slice() = terminalAugQ;
                
        dmat tempK(K_old.n_rows, K_old.n_cols);
        dmat tempL(L_old.n_rows, L_old.n_cols);
        dmat tempM(M_old.n_rows, M_old.n_cols);
                
        // backward pass to compute the gains
        for (auto k= horizonLength-1; k >= 0; --k) {
            // evaluate the epsilon operator
            costate currentCostateEpsilon = evaluateEpsilonOperator(currentCostate, transitionMatrixCa);
            
            computeGains(Xu.at(k), Xl.at(k), currentCostateEpsilon, modeProbsCa.colptr(k), transitionMatrixCa, augA, augB, augC, 
                augV, Se.slice(k), JRJ, tempK, tempL, tempM);            
        
            const costate tempCostate = computeCostate(currentCostateEpsilon, augA, augB, tempK, tempL, tempM, augV, augW, Se.slice(k), augC,
                transitionMatrixCa, augQ, JRJ);
         
            double costToGo = dot(tempCostate.omega, modeProbsCa.col(k));
            for (auto i = 0; i < numCaModes; ++i) {
                costToGo += trace(tempCostate.Pu.slice(i) * (Xu.at(k).slice(i) + Xl.at(k).slice(i)) 
                    + tempCostate.Pl.slice(i) * Xu.at(k).slice(i));
            }            
            if (costToGo < oldCostToGo.at(k)) {
                oldCostToGo(k) = costToGo;
                K_new.slice(k) = tempK;
                L_new.slice(k) = tempL;
                M_new.slice(k) = tempM;
                currentCostate = tempCostate;
            } else {
                // recompute the costate with the old gain
                currentCostate = computeCostate(currentCostateEpsilon, augA, augB, K_new.slice(k), L_new.slice(k), M_new.slice(k), 
                    augV, augW, Se.slice(k), augC, transitionMatrixCa, augQ, JRJ);                            
            }            
        }
        if (abs(oldCost - oldCostToGo(0)) < CONVERGENCE_DIFF) {            
            break;
        } else {
            oldCost = oldCostToGo(0);
        }        
    }
    
    // compute the new moments
    outputArrays[3] = armaCreateMxMatrix(augStateCov_old.n_rows, augStateCov_old.n_cols, numCaModes); // augStateCov
    outputArrays[4] = armaCreateMxMatrix(augStateCov_old.n_rows, augStateCov_old.n_cols, numCaModes); // augStateSecondMoment
            
    dcube augStateCov_new(mxGetPr(outputArrays[3]), augStateCov_old.n_rows, augStateCov_old.n_cols, numCaModes, 
            false, true);
    dcube augStateSecondMoment_new(mxGetPr(outputArrays[4]), augStateCov_old.n_rows, augStateCov_old.n_cols, numCaModes, 
            false, true);
    
    const dmat KS = K_new.slice(0) * Se.slice(0);
    const dmat KSC = KS * augC;
    const dmat MKSC = M_new.slice(0) + KSC;
    const dmat E_tilde = symmatu(KS * augV * KS.t());
    const double* modeCol = modeProbsCa.colptr(0);

    dcube firstPart(E_tilde.n_rows, E_tilde.n_cols, numCaModes);
    dcube secondPart(E_tilde.n_rows, E_tilde.n_cols, numCaModes);
    //#pragma omp parallel for shared(AKSC, AKSCMBL)
    for (auto i = 0; i < numCaModes; ++i) {
        const dmat AKSC = augA.slice(i) - KSC;
        const dmat D_tilde = AKSC  + augB.slice(i) * L_new.slice(0) - M_new.slice(0); // D_tilde in the paper
        firstPart.slice(i) = symmatu(AKSC * augStateCov_old.slice(i) * AKSC.t() 
            + D_tilde * augStateSecondMoment_old.slice(i) * D_tilde.t())
            + modeCol[i] * (E_tilde + augW);
        secondPart.slice(i) = symmatu(MKSC * augStateSecondMoment_old.slice(i) * MKSC.t()
            + KSC * augStateCov_old.slice(i) * KSC.t())
            + modeCol[i] * E_tilde;
    }

    //#pragma omp parallel for shared(modeCol)
    for (auto j=0; j < numCaModes; ++j) { // new mode            
        const double* col = transitionMatrixCa.colptr(j);
        for (auto i=0; i < numCaModes; ++i) { // old mode
            augStateCov_new.slice(j) += col[i] * firstPart.slice(i);
            augStateSecondMoment_new.slice(j) += col[i] * secondPart.slice(i); 
        }
    }                       
    
    outputArrays[5] = mxCreateDoubleScalar(static_cast<double>(counter)); //numIterations
}

