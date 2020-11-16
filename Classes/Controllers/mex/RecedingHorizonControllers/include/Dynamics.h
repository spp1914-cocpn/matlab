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

#ifndef _RECEDING_HORIZON_UDPLIKE_CONTROLLER_DYNAMICS_H
#define _RECEDING_HORIZON_UDPLIKE_CONTROLLER_DYNAMICS_H

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
        dmat& modeProbs;
        dcube& Se; 
        
    public:
        Dynamics(dcube& augA, dcube& augB, sp_dmat& augC, dmat& augW, dmat& augV,
                dmat& transitionMatrixCa, dmat& modeProbsCa, dcube& Se) 
         : A(augA), B(augB), C(augC), W(augW), V(augV), transitionMatrix(transitionMatrixCa),
            modeProbs(modeProbsCa), Se(Se) {
         }  
         
         ~Dynamics() {}
         
         Dynamics(const Dynamics& other) 
            : A(other.A), B(other.B), C(other.C), W(other.W), V(other.V), transitionMatrix(other.transitionMatrix),
            modeProbs(other.modeProbs), Se(other.Se) {}
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
            this->modeProbs = other.modeProbs;
            this->Se = other.Se;

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
         
         const int getNumModes() const {
             return this->transitionMatrix.n_rows;
         }
         
         const int getHorizonLength() const {
             return this->Se.n_slices;
         }
    };
}
#endif
