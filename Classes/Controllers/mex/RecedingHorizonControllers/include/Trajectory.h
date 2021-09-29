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

#ifndef _RECEDING_HORIZON_UDPLIKE_CONTROLLER_TRAJECTORY_H
#define _RECEDING_HORIZON_UDPLIKE_CONTROLLER_TRAJECTORY_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <armadillo>
#include "Dynamics.h"

using namespace arma;

namespace RecedingHorizonUdpLikeController {
    
    class Trajectory {
        
        private:
            const Dynamics& dynamics;
            const int horizonLength;
                        
            field<dcube> Xu;
            field<dcube> Xl;         
                    
        private:
            /* Non-copyable*/
            Trajectory(const Trajectory& other) = delete;
            Trajectory(const Trajectory&& other) = delete;
            Trajectory operator=(const Trajectory& other) = delete;
            
        public:
            Trajectory(const Dynamics& dyn, const dcube& augStateCov, const dcube& augStateSecondMoment) 
                : dynamics(dyn), horizonLength(dyn.getHorizonLength()) {
                
                this->Xu = field<dcube>(horizonLength + 1);
                this->Xl = field<dcube>(horizonLength + 1);
                
                this->Xu.for_each([&augStateCov](dcube& c) {c.copy_size(augStateCov);});
                this->Xl.for_each([&augStateSecondMoment](dcube& c) {c.copy_size(augStateSecondMoment);});
                this->Xu.at(0) = augStateCov;
                this->Xl.at(0) = augStateSecondMoment;                
            }
                
            ~Trajectory() {}
            
        
            const dcube& getXu(const int k) const {
                return this->Xu.at(k);
            }
            
            const dcube& getXl(const int k) const {
                return this->Xl.at(k);
            }
            
            double computeCostToGoForCostate(const int stage, const Costate& forCostate) const {
                double costToGo = dot(forCostate.getOmega(), this->dynamics.getModeProbsAtCol(stage));
                
                for (auto i = 0; i < this->dynamics.getNumModes(); ++i) {
                    costToGo += trace(forCostate.getPu().slice(i) * (this->getXu(stage).slice(i) + this->getXl(stage).slice(i)) 
                        + forCostate.getPl().slice(i) * this->getXu(stage).slice(i));
                }
                return costToGo;
            }
            
            void predictHorizon(const dcube& K, const dcube& L) {
                
                const int numModes = this->dynamics.getNumModes();
                
                // the first element is already filled
                 for (auto k= 0; k < this->horizonLength; ++k) {
                    const dmat KS = K.slice(k) * this->dynamics.getSeAt(k);
                    const dmat KSC = KS * this->dynamics.getAugC();
                    const dmat E_tilde = symmatu(KS * this->dynamics.getAugV() * KS.t()); // E_tilde in the paper
                    const dmat noisePartFull = E_tilde + this->dynamics.getAugW();
                    
                    const dmat AexpBexpL = this->dynamics.getAexpAt(k) + this->dynamics.getBexpAt(k) * L.slice(k); //Ahat_k+Bhat_k*L_k

                    const dcube currXu = this->Xu.at(k);
                    const dcube currXl = this->Xl.at(k);
                    const double* modeCol = this->dynamics.getModeProbsAt(k);
                    dcube firstPart(E_tilde.n_rows, E_tilde.n_cols, numModes);
                    dcube secondPart(E_tilde.n_rows, E_tilde.n_cols, numModes);
                    
                    #ifdef _OPENMP
                    #pragma omp parallel for shared(firstPart, secondPart)
                    #endif
                    for (auto i = 0; i < numModes; ++i) {
                        const dmat AKSC = this->dynamics.getAugA(i) - KSC;
                        const dmat diffPart = this->dynamics.getAugA(i) + this->dynamics.getAugB(i) * L.slice(k) - AexpBexpL;
                        //const dmat diffPart = this->dynamics.getAugA(i) - this->dynamics.getAexpAt(k) + (this->dynamics.getAugB(i) - this->dynamics.getBexpAt(k)) * L.slice(k);
                                               
                        firstPart.slice(i) = symmatu(AKSC * currXu.slice(i) * AKSC.t() + diffPart * currXl.slice(i) * diffPart.t())
                                        + modeCol[i] * noisePartFull;
                        secondPart.slice(i) = symmatu(AexpBexpL * currXl.slice(i) * AexpBexpL.t() + KSC * currXu.slice(i) * KSC.t())
                                        + modeCol[i] * E_tilde;
                    }
                    this->Xu.at(k+1).zeros();
                    this->Xl.at(k+1).zeros();
        
                    #ifdef _OPENMP        
                    #pragma omp parallel for
                    #endif
                    for (auto j=0; j < numModes; ++j) { // new mode            
                        const double* col = this->dynamics.getTransitionMatrixCol(j);
                        for (auto i=0; i < numModes; ++i) { // old mode
                            this->Xu.at(k+1).slice(j) += col[i] * firstPart.slice(i);
                            this->Xl.at(k+1).slice(j) += col[i] * secondPart.slice(i);
                        }            
                    }
                }                
            }
    };
}
#endif
