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
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace arma;

// implements a vectorized element-wise mode operation: res[i] = mod(a[i], m)
inline ucolvec arma_mod(const ucolvec& a, const int m) {
    return a - floor(a / m) * m;
}

// computes the stationary distribution of the augmented delay chain, based on the stationary distribution of the delay chain
const dvec computeStationaryDistribution(const dmat& delayTransitionMatrix, const umat& statesBaseN, 
        const uword numAggregations, const uword numAggregatedStates) {
    // shortcut in case numAggregations is N=2:
    // compute stationary distribution of aggT by using stationary distribution of delayTransitionMatrix
    // by employing Theorem 6.5.2 in Kemeny and Snell
    // could be extended to the general case as well
    dvec statDist;
    cx_vec eigvals;
    cx_mat reigvecs; // we need the left eigenvectors, which are the right eigenvectors of the transposed matrix
    eig_gen(eigvals, reigvecs, delayTransitionMatrix.t());
    // find the eigenvector corresponding to the largest eigenvalue 1 (Perron eigenvalue) and normalize (entries sum to 1)
    const dcolvec eigVec = conv_to<dcolvec>::from(reigvecs.col(index_max(eigvals))) / sum(conv_to<dcolvec>::from(reigvecs.col(index_max(eigvals))));
    if (numAggregations == 2) {
        statDist = (delayTransitionMatrix.each_col() % eigVec).as_row();        
    } else {        
        statDist.set_size(numAggregatedStates, 1);
        //dcolvec statDist2 = zeros<dcolvec>(numAggregatedStates);
        // [N N-1 N-1 N-2 N-2 ... 1 1 0]
        //const urowvec rowIdx = regspace<urowvec>(numAggregations-1, -1, 1);
        //const urowvec colIdx = regspace<urowvec>(numAggregations-2, -1, 0);        
        #ifdef _OPENMP    
        #pragma omp parallel for shared(statDist)
        #endif
        for (uword i = 0; i < numAggregatedStates; ++i) {
            const urowvec currStateBaseN = statesBaseN.row(i);
            double prod = eigVec(currStateBaseN(numAggregations-1)); // last element in row
            // straightforward implementation
            for (uword j= numAggregations-1; j > 0; --j) {
                prod *= delayTransitionMatrix(currStateBaseN(j), currStateBaseN(j-1)); 
            }
            statDist(i) = prod;
            /*const umat subIdx = join_rows(currStateBaseN.elem(rowIdx), currStateBaseN.elem(colIdx)).t();
            const double prod = eigVec(currStateBaseN(numAggregations-1)) 
             * arma::prod(delayTransitionMatrix.elem(sub2ind(size(delayTransitionMatrix), subIdx)));
            statDist(i) = prod;*/
        }        
    }
    return statDist;
}

// fills the distribution matrix U and the collection matrix V for the lumped dynamics P = U*T_agg*V
// the entries of the distribution matrix are 1/d_j with d_j the number of elements in cluster j
void fillMatricesVU(umat& V, dmat& U,  const umat& statesBaseN, const uword numDelays) {
    
    auto numCaModes = V.n_cols;
    auto numAggregatedStates = V.n_rows;
    // first cluster is easy: tau_k = 0, which is regularly spaced, starting at zero
    const uvec clusterIdx = regspace<ucolvec>(0, numDelays, numAggregatedStates - 1);    
    U.submat(uvec({0}), clusterIdx).fill(1.0 / (numAggregatedStates/numDelays));
    V.submat(clusterIdx, uvec({0})).ones();
    // i-th cluster: tau_k > 0, tau_{k-1} > 1, ... tau_{k-(i-1)} > i-1, tau_{k-i} <= i
    // last cluster: tau_k > 0, tau_{k-1} > 1, ... tau_{k-(i-1)} > i-1, tau_{k-i} > i
    #ifdef _OPENMP    
    #pragma omp parallel for shared(U,V)
    #endif
    for (uword i = 1; i < numCaModes; ++i) {        
        ucolvec allFoundIdx;
        ucolvec foundIdx;
        // check for the last cluster, has to be treated differently
        if (i == numCaModes - 1) {
            foundIdx = find(all(statesBaseN.head_cols(i-1) > repmat(regspace<urowvec>(0, 1, i-2), numAggregatedStates, 1), 1));
            allFoundIdx = intersect(foundIdx, find(statesBaseN.unsafe_col(i-1) > i-1));  // tau_{k-i} > i
        }
        else {
            // tau_k > 0, tau_{k-1} > 1, ... tau_{k-(i-1)} > i-1
            foundIdx = find(all(statesBaseN.head_cols(i) > repmat(regspace<urowvec>(0, 1, i-1), numAggregatedStates, 1), 1));
            allFoundIdx = intersect(foundIdx, find(statesBaseN.unsafe_col(i) <= i));  // tau_{k-i} <= i    
        }
        U.submat(uvec({i}), allFoundIdx).fill(1.0 / allFoundIdx.n_elem);
        V.submat(allFoundIdx, uvec({i})).ones();
    }
}

// fills the distribution matrix U and the collection matrix V for the lumped dynamics P = U*T_agg*V
// the distribution matrix is specified by the stationary distribution of aggT
void fillMatricesVU(umat& V, dmat& U,  const umat& statesBaseN, const uword numDelays, const dvec& statDist) {
    
    auto numCaModes = V.n_cols;
    auto numAggregatedStates = V.n_rows;
    // first cluster is easy: tau_k = 0, which is regularly spaced, starting at zero
    const uvec clusterIdx = regspace<ucolvec>(0, numDelays, numAggregatedStates - 1);    
    //U.submat(uvec({0}), clusterIdx).fill(1.0 / (numAggregatedStates/numDelays));
    U.submat(uvec({0}), clusterIdx) = conv_to<drowvec>::from(statDist.elem(clusterIdx) / sum(statDist.elem(clusterIdx)));
    V.submat(clusterIdx, uvec({0})).ones();
    // i-th cluster: tau_k > 0, tau_{k-1} > 1, ... tau_{k-(i-1)} > i-1, tau_{k-i} <= i
    // last cluster: tau_k > 0, tau_{k-1} > 1, ... tau_{k-(i-1)} > i-1, tau_{k-i} > i
    #ifdef _OPENMP    
    #pragma omp parallel for shared(U,V)
    #endif
    for (uword i = 1; i < numCaModes; ++i) {        
        ucolvec allFoundIdx;
        ucolvec foundIdx;
        // check for the last cluster, has to be treated differently
        if (i == numCaModes - 1) {
            foundIdx = find(all(statesBaseN.head_cols(i-1) > repmat(regspace<urowvec>(0, 1, i-2), numAggregatedStates, 1), 1));
            allFoundIdx = intersect(foundIdx, find(statesBaseN.unsafe_col(i-1) > i-1));  // tau_{k-i} > i
        }
        else {
            // tau_k > 0, tau_{k-1} > 1, ... tau_{k-(i-1)} > i-1
            foundIdx = find(all(statesBaseN.head_cols(i) > repmat(regspace<urowvec>(0, 1, i-1), numAggregatedStates, 1), 1));
            allFoundIdx = intersect(foundIdx, find(statesBaseN.unsafe_col(i) <= i));  // tau_{k-i} <= i    
        }
        //U.submat(uvec({i}), allFoundIdx).fill(1.0 / allFoundIdx.n_elem);
        U.submat(uvec({i}), allFoundIdx) = conv_to<drowvec>::from(statDist.elem(allFoundIdx) / sum(statDist.elem(allFoundIdx)));
        V.submat(allFoundIdx, uvec({i})).ones();      
    }
}

void mexFunction(int numOutputs, mxArray* outputArrays[], 
        int numInputs, const mxArray* inputArrays[]) {
    
    // mandatory arguments: transition matrix of the delay MC tau_k and the desired sequence length
    const dmat delayTransitionMatrix = armaGetPr(inputArrays[0]);
    const uword sequenceLength = armaGetScalar<uword>(inputArrays[1]);
    
    const uword numDelays = delayTransitionMatrix.n_rows; // last state also comprises packet losses
    const uword numCaModes = sequenceLength + 1; // number of clusters/lumps
    
    // optional argument, defaults to true
    const bool useStationaryDist = (numInputs > 2 && !mxIsEmpty(inputArrays[2])) ? armaGetScalar<bool>(inputArrays[2]) : true;
    
    if (!delayTransitionMatrix.is_square() || numDelays < numCaModes) {
        // something went wrong
        mexErrMsgIdAndTxt("Utility:mex_CalculateTransitionMatrixCorrDelays:InvalidTransitionMatrix", 
                "** Transition matrix of the delay Markov chain must be at least %d-by-%d **", numCaModes, numCaModes);
        return;
    }    
    
    // we return the transition matrix of the lumped chain
    outputArrays[0] = mxCreateUninitNumericMatrix(static_cast<size_t>(numCaModes), static_cast<size_t>(numCaModes), 
            mxDOUBLE_CLASS, mxREAL);    
    dmat transitionMatrixCa(mxGetPr(outputArrays[0]), numCaModes, numCaModes, false, true);
    
    const uword numAggregations = sequenceLength;
    const uword numAggregatedStates = pow(numDelays, numAggregations);
    // the aggregated chain is not regular
    const ucolvec aggregatedStatesDecimal = regspace<ucolvec>(0, numAggregatedStates - 1);
    umat statesBaseN = zeros<umat>(numAggregatedStates, numAggregations);
    
    #ifdef _OPENMP    
    #pragma omp parallel for shared(statesBaseN)
    #endif
    for (uword i=0; i < numAggregations; ++i) {
        const ucolvec a = floor(aggregatedStatesDecimal / pow(numDelays, i));
        statesBaseN.unsafe_col(i) = arma_mod(a, numDelays); // implements vectorized mod(a, numDelays)
    }
    
    // matrix has repetitive cations =structure, so only compute first fraction
    const uword numRows = numAggregatedStates / numDelays;
    // use a sparse matrix for the transitions of the aggregated chain
    const int numNonzeroEntries = std::pow(numDelays, numAggregations); // maximum number of nonzero transition probs of the aggregated chain
    // compute the locations of the sparse entries
    const umat sparseLocations = join_vert(kron(regspace<urowvec>(0, numRows - 1), ones<urowvec>(numDelays)), regspace<urowvec>(0, numNonzeroEntries - 1));
    // values of the aggregated matrix are the repeated entries of the original one
    const dcolvec values = repmat(delayTransitionMatrix.as_row(), 1, numNonzeroEntries / delayTransitionMatrix.n_elem).as_col();
    
    const sp_dmat aggT = repmat(sp_dmat(sparseLocations, values, numRows, numAggregatedStates), numDelays, 1);
              
    // compute the indices of the clusters: states of aggregated chain are lumped into clusters
    // -> mode of augmented system
    // construct matrices U, V required for clustering into modes (theta_k) (MC lumping)
    dmat U = zeros<dmat>(numCaModes, numAggregatedStates);
    umat V = zeros<umat>(numAggregatedStates, numCaModes); // binary matrix
    
    if (useStationaryDist) {
        const dvec statDist = computeStationaryDistribution(delayTransitionMatrix, statesBaseN, numAggregations, numAggregatedStates);
        fillMatricesVU(V, U,  statesBaseN, numDelays, statDist);
    } else {
        fillMatricesVU(V, U,  statesBaseN, numDelays);
    }
    
    // in general, aggT is not exactly lumpable, so this is just an approximation
    // but this choice of U ensures that aggT (being primitive/irreducible) and transitionMatrixCa have the "same"
    // stationary distribution (w.r.t the clusters)
    // meaning that the dynamics of the original chain, constrained to the clusters (i.e. theta_k), has the
    // same asymptotic behavior as the lumped dynamics
    transitionMatrixCa = U * aggT * V;     
}
