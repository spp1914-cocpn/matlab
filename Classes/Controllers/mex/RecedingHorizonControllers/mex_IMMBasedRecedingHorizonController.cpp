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
#pragma omp declare reduction(+ : dmat : omp_out += omp_in) initializer(omp_priv = zeros<dmat>(size(omp_orig)))
#pragma omp declare reduction(+ : dcolvec : omp_out += omp_in) initializer(omp_priv = omp_orig)
#endif

using namespace std;
using namespace arma;

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // required for computation of input sequence: posterior means and mode probs, eta, augA, augB, costate, transition matrix
    // and augR
    const dcube augA = armaGetCubePr(inputArrays[0]);
    const dcube augB = armaGetCubePr(inputArrays[1]);
    const dmat transitionMatrix = armaGetPr(inputArrays[2]);
    const dmat augR = armaGetPr(inputArrays[3]); // only for first mode
    const dmat stateMeans = armaGetPr(inputArrays[4]);
    const dcolvec modeProbs = armaGetPr(inputArrays[5]);
    const dcolvec eta = armaGetPr(inputArrays[6]);
    const dcube costate = armaGetCubePr(inputArrays[7]);
    
    const uword numModes = augA.n_slices;
    const uword dimAugState = augA.n_rows;
    const uword dimInputSequence = augB.n_cols;
    const uword dimPlantState = stateMeans.n_rows;
    const uword dimEta = dimAugState - dimPlantState;
    
    outputArrays[0] = mxCreateUninitNumericMatrix(static_cast<size_t>(dimInputSequence), 1, mxDOUBLE_CLASS, mxREAL);
    dcolvec inputSeq(mxGetPr(outputArrays[0]), dimInputSequence, false, true);
    
    dmat L = modeProbs.at(0) * augR;
    dcolvec part2 = zeros<dcolvec>(dimInputSequence);
    const dmat A = augA(span(0, dimPlantState - 1), span(0, dimPlantState - 1), span(0));
    const dmat F = augA(span(dimAugState - dimEta, dimAugState - 1), span(dimAugState - dimEta, dimAugState - 1), span(0));
    const dcolvec etaPred = F * eta; // F*eta, same for all modes
    
    #ifdef _OPENMP    
    #pragma omp parallel for reduction(+:L, part2)
    #endif
    for (uword i= 0; i < numModes; ++i) {
        const dmat BH = augA(span(0, dimPlantState - 1), span(dimAugState - dimEta, dimAugState - 1), span(i));
        const dcolvec statePart = join_cols(A * stateMeans.col(i) + BH * eta, etaPred);
        //const dcolvec mean = join_cols(stateMeans.col(i), eta);
        mat::const_row_iterator row = transitionMatrix.begin_row(i);        
        dmat weightedCostate = zeros<dmat>(dimAugState, dimAugState);
        for (uword r = 0; r < numModes; ++r, ++row) {            
            weightedCostate += (*row) * costate.slice(r);
        }
        part2 += modeProbs.at(i) * augB.slice(i).t() * weightedCostate * statePart;
        L += modeProbs.at(i) * augB.slice(i).t() * weightedCostate * augB.slice(i);
    }
    inputSeq = solve(symmatu(L), -part2);
}
