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

#ifndef _RECEDING_HORIZON_UDPLIKE_CONTROLLER_COSTATE_H
#define _RECEDING_HORIZON_UDPLIKE_CONTROLLER_COSTATE_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include "Dynamics.h"
#include <armadillo>

using namespace arma;

namespace RecedingHorizonUdpLikeController {
    
    class Costate {
        
        private:
            dcube Pu;
            dcube Pl;
            dcolvec omega;
            
            int stage;
            Dynamics& dynamics;
        
        private:
            Costate(const dcube& Pu, const dcube& Pl, const dcolvec& omega, int k, Dynamics& dyn) 
                : Pu(Pu), Pl(Pl), omega(omega), stage(k), dynamics(dyn) {}                        
                
        public:
            Costate(const Costate& other) 
                : Pu(other.Pu), Pl(other.Pl), omega(other.omega), stage(other.stage), dynamics(other.dynamics) {}
            
            Costate& operator=(const Costate& other) {
                if (this == &other) {
                    return *this;
                }    
                this->Pu = other.Pu;
                this->Pl = other.Pl;
                this->omega = other.omega;
                this->dynamics = other.dynamics;
                this->stage = other.stage;                
  
                return *this;
            }           
            
            ~Costate() {}
            
            const dcube& getPu() const {
                return this->Pu;
            }
            
            const dcube& getPl() const {
                return this->Pl;
            }
            
            const dcolvec& getOmega() const {
                return this->omega;
            }
            
            Costate evaluateEpsilonOperator() const {
                Costate result(zeros<dcube>(size(this->Pu)), zeros<dcube>(size(this->Pl)), 
                        this->dynamics.getTransitionMatrix() * this->omega, this->stage, this->dynamics);
                
                #ifdef _OPENMP        
                #pragma omp parallel for shared (result)
                #endif
                for (auto i=0; i < this->dynamics.getNumModes(); ++i) {
                    mat::const_row_iterator row = this->dynamics.getTransitionMatrix().begin_row(i);
                    for (auto j=0; j < this->dynamics.getNumModes(); ++j, ++row) {
                        result.Pu.slice(i) += (*row) * this->Pu.slice(j);
                        result.Pl.slice(i) += (*row) * this->Pl.slice(j);            
                    }  
                }
                return result;
            }           
            
            Costate computeNew(const dmat& K, const dmat& L, const dcube& augQ, const dmat& JRJ) const {
                
                Costate result(augQ, zeros<dcube>(size(augQ)), this->dynamics.getTransitionMatrix() * this->omega,
                        this->stage - 1, this->dynamics);

                result.Pu.slice(0) += symmatu(L.t() * JRJ * L);
                result.Pl.slice(0) += symmatu(L.t() * JRJ * L);
                                               
                const dmat KS = K * this->dynamics.getSeAt(result.stage);
                const dmat KSC = KS * this->dynamics.getAugC();
                const dmat E_tilde = symmatu(KS * this->dynamics.getAugV() * KS.t()); // E_tilde in the paper            
                const dmat AexpBexpL = this->dynamics.getAexpAt(result.stage) +this->dynamics.getBexpAt(result.stage)  * L;
                const dmat AexpBexpLKSC = AexpBexpL - KSC; //Aexp+Bexp*L-K*S*C
                
                #ifdef _OPENMP
                #pragma omp parallel for shared(result)
                #endif
                for (int i = 0; i < this->dynamics.getNumModes(); ++i) {
                    const dmat currentPlEpsilon = this->Pl.slice(i);
                    const dmat currentPuEpsilon = this->Pu.slice(i);
                    const dmat BL = this->dynamics.getAugB(i) * L;
                    const dmat ABL = this->dynamics.getAugA(i) + BL; // U_tilde in the paper
                    const dmat diffPart = ABL - AexpBexpL;   
                    const dmat O_tilde = AexpBexpLKSC - BL; //Aexp+Bexp*L-K*S*C-B(i)*L

                    result.omega(i) += trace(currentPlEpsilon * E_tilde + (currentPlEpsilon + currentPuEpsilon) * this->dynamics.getAugW());
                    result.Pu.slice(i) += symmatu(ABL.t() * currentPuEpsilon * ABL + diffPart.t() * currentPlEpsilon * diffPart);
                    result.Pl.slice(i) += symmatu(BL.t() * currentPuEpsilon * BL + O_tilde.t() * currentPlEpsilon * O_tilde);        
                }
                
                return result;
            }            
   
            static Costate createTerminalCostate(Dynamics& dyn, const dmat& terminalAugQ) {
                dcube terminalPu = dcube(terminalAugQ.n_rows, terminalAugQ.n_cols, dyn.getNumModes());
                terminalPu.each_slice() = terminalAugQ;

                Costate terminalCostate(terminalPu, zeros<dcube>(size(terminalPu)), zeros<dcolvec>(dyn.getNumModes()),
                        dyn.getHorizonLength(), dyn);

                return terminalCostate;
            }
    };
}
#endif
