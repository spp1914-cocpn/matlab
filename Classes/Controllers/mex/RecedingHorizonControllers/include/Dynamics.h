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

#ifndef _RECEDING_HORIZON_UDPLIKE_CONTROLLER_DYNAMICS_H
#define _RECEDING_HORIZON_UDPLIKE_CONTROLLER_DYNAMICS_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <armadillo>

using namespace arma;
using namespace std;

namespace RecedingHorizonUdpLikeController {
    
    class Dynamics {
        
    private:
        dcube& A;
        dcube& B;        
        sp_dmat& C;
        dmat& W;
        dmat& V;

        dmat& transitionMatrix;                        
        dcube& Se;  
           
        dmat modeProbs;
        dcube Aexp;
        dcube Bexp;
        
    public:
        Dynamics(dcube& augA, dcube& augB, sp_dmat& augC, dmat& augW, dmat& augV,
                dmat& transitionMatrixCa, dcolvec& modeProbsCa, dcube& Se) 
         : A(augA), B(augB), C(augC), W(augW), V(augV), transitionMatrix(transitionMatrixCa), Se(Se){
                        
            // compute Aexp and Bexp for the horizon
            // also the mode probs
            this->Aexp = zeros<dcube>(this->A.n_rows, this->A.n_rows, this->getHorizonLength());
            this->Bexp = zeros<dcube>(this->A.n_rows, this->B.n_cols, this->getHorizonLength());           
            
            this->modeProbs = join_rows(modeProbsCa, zeros<dmat>(this->getNumModes(), this->getHorizonLength()));
            
            // compute Aexp and Bexp for the horizon            
            for (auto k=0; k < this->getHorizonLength(); ++k) {                
                const dcolvec probs = this->getModeProbsAtCol(k);        
                for (auto i=0; i < this->getNumModes(); ++i) {
                    this->Aexp.slice(k) += probs[i] * this->A.slice(i);
                    this->Bexp.slice(k) += probs[i] * this->B.slice(i);
                }
                // predict the mode (theta) probabilities
                this->modeProbs.col(k+1) = this->transitionMatrix.t() * probs;
            }
         }  
         
         ~Dynamics() {}
         
         Dynamics(const Dynamics& other) 
            : A(other.A), B(other.B), C(other.C), W(other.W), V(other.V), transitionMatrix(other.transitionMatrix), 
                    Se(other.Se) {        
            
            this->modeProbs = other.modeProbs;
            this->Aexp = other.Aexp;
            this->Bexp = other.Bexp;            
        }
//          
         Dynamics& operator=(const Dynamics& other) {
            if (this == &other) {
                return *this;
            }
            this->A = other.A;
            this->B = other.B;
            this->C = other.C;
            this->W = other.W;
            this->V = other.V;
            this->transitionMatrix = other.transitionMatrix;            
            this->Se = other.Se;
            
            this->modeProbs = other.modeProbs;
            this->Aexp = other.Aexp;
            this->Bexp = other.Bexp;

            return *this;
         }
         
         const dmat& getAugA(int slice) const {
             return this->A.slice(slice);
         }
         
         const dmat& getAugB(int slice) const {
             return this->B.slice(slice);
         }
         
         const sp_dmat& getAugC() const {
            return this->C;
         }
         
         const dmat& getAugW() const {
             return this->W;
         }
         
         const dmat& getAugV() const {
             return this->V;
         }
         
         const dmat& getTransitionMatrix() const {
             return this->transitionMatrix;
         }
         
         const double* getTransitionMatrixCol(int k) const {
             return this->transitionMatrix.colptr(k);
         }
         
         const dmat& getModeProbs() const {
             return this->modeProbs;
         }
         
         const double* getModeProbsAt(int k) const {
             return this->modeProbs.colptr(k);
         }
         
         const dcolvec getModeProbsAtCol(int k) const {
             return this->modeProbs.col(k);
         }
         
         const dcube& getSe() const {
             return this->Se;
         }
         
         const dmat& getSeAt(int k) const {
             return this->Se.slice(k);
         }
         
         const dmat& getAexpAt(int k) const {
             return this->Aexp.slice(k);
         }
         
         const dmat& getBexpAt(int k) const {
             return this->Bexp.slice(k);
         }
         
         const int getNumModes() const {
             return this->transitionMatrix.n_rows;
         }
         
         const int getHorizonLength() const {
             return this->Se.n_slices;
         }
    };
}
#endif
