classdef BallBeamSystem < NonlinearPlant
    % Implementation of a ball rolling freely on a rotating beam with  
    % nonlinear dynamics according to x_dot = f(x,u,w)
    % x=[q;q_dot;i] is the 5-dimensional state vector consisting of
    % the generalized coordinates q =[r; theta], consisting of the displacement of the ball (r, in meters) along the beam 
    % and the angular rotation of the beam (theta, in radians), and the
    % current flowing through the servo-motor (i, in amperes) that generates the torque required to rotate the beam up and down. 
    % u is the control input (voltage, in volts) applied to the servo-motor to produce the torque,
    % and w is an optional process noise, modeling disturbance forces (in
    % newtons and volts, respectively) acting on the ball's movement (i.e,
    % in direction of r) and the control input.
    % Note that the r and theta wer chosen such that r=0 represents the center of the beam, 
    % r>0 indicates a displacement to the right, and theta > 0 means a counter-clockwise beam rotation. 
    % Consequently, the unstable equilibrium is given by x=[0;0;0;0;0].
        
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2018-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
    %
    %                        Institute for Anthropomatics and Robotics
    %                        Chair for Intelligent Sensor-Actuator-Systems (ISAS)
    %                        Karlsruhe Institute of Technology (KIT), Germany
    %
    %                        https://isas.iar.kit.edu
    %
    %    This program is free software: you can redistribute it and/or modify
    %    it under the terms of the GNU General Public License as published by
    %    the Free Software Foundation, either version 3 of the License, or
    %    (at your option) any later version.
    %
    %    This program is distributed in the hope that it will be useful,
    %    but WITHOUT ANY WARRANTY; without even the implied warranty of
    %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %    GNU General Public License for more details.
    %
    %    You should have received a copy of the GNU General Public License
    %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (Constant, Access = private)
        graviationalAcceleration = 9.81; % in N/kg
        defaultInternalSamplingRate = 0.001; % 1 kHz
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        massBeam; % M, mass of beam, located at d meters from pivot, in kg
        lengthBeam; % L, in m
        distancePivotPlaneBallContact; % d, distance from pivot to plane of ball contact on beam in m, >= distance from pivot to beam center of mass
        distancePivotBeamCenter; % D, distance from pivot to beam center of mass, in m
        
        massBall; % m, in kg
        radiusBall; % R_0 radius of the ball, in m
        distanceBallAxisBeam; % R_1, distance between ball's axis of rotation and its point of contact with beam, in m, <= ball radius     
        
        inertiaBeam; % J_beam, beam's moment of inertia, including beam mass offset and all rotating components, in kgm^2
        
        frictionBeamBall; % C_1, viscous friction coefficient between ball and beam, in Ns/m
        frictionMotorBeam; % C_2, viscous friction coefficient/damping constant of servo-motor, in Nm/(rad/s))
        % Nm/(rad/s) is Nms
        % typical value: 0.4 * 10^-6
        
        electricalResistance; % R_m, servo-motor's electrical resistance, in ohms
        % typical value: 4-5 ohms
        electricalInductance; % L_m, servo-motor's electrical inductance, in henrys
        % typically very small, less than 1 H, e.g. 0.11mH
        motorTorqueConstant; % K_m, servo-motor's electromotive torque constant, equals its back emf constant, in Nm/A        
        % typical value: 0.01 Nm/A
        
        % we measure theta and r directly
        % or alternatively, i and r
        % typical parameter values: m = 0.12 kg, M = 0.5kg, R_0 = 0.015m
        % ball radius
        % L = 0.8m
    end
    
    properties (SetAccess = private, GetAccess = public)
        %
        % variance of disturbance force (in newtons) driving the movement of the ball
        % default zero
        varDisturbanceForceBall (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        % variance of disturbance (in volts) driving the motor (entering additively at the actuator, actuation noise)
        % default zero
        varDisturbanceForceActuator (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
        %        
    end
    
    properties (Access = public)
        samplingInterval(1,1) double {mustBePositive, mustBeFinite} = BallBeamSystem.defaultInternalSamplingRate;
        
        % we need noise terms for the continuous-time linearization of the model
        % can be different from the ones used to simulate the nonlinear dynamics
        % the controllers use discrete-time linearizations, the noise in
        % these models is dependent on the sampling rate
        % the terms are assumed to be independent of each other
        %
        varDisturbanceForceBallContLin (1,1) double {mustBeNonnegative, mustBeFinite} = 0;        
        varDisturbanceForceActuatorContLin (1,1) double {mustBeNonnegative, mustBeFinite} = 0;
    end
    
    methods (Access = public)
        %% BallBeamSystem
        function this = BallBeamSystem(massBeam, lengthBeam, inertiaBeam, massBall, radiusBall, ...
                distancePivotPlaneBallContact, distancePivotBeamCenter, distanceBallAxisBeam, ...
                frictionBeamBall, frictionMotorBeam, electricalResistance, electricalInductance, motorTorqueConstant, samplingInterval)
            % Class constructor.
            %
            % Parameters:
            %   >> massBeam (Positive Scalar)
            %      A positive scalar denoting the mass (in kg) of the beam, M.
            %
            %   >> lengthBeam (Positive scalar)
            %      A positive scalar denoting the length (in m) of the beam, L.
            %
            %   >> inertiaBeam (Positive scalar)
            %      A positive scalar denoting the moment of inertia (in kgm^2) of the beam, J_beam, 
            %      including all rotational components and including the offset of the beam mass.
            %
            %   >> massBall (Positive scalar)
            %      A positive scalar denoting the mass (in kg) of the ball, m.          
            %
            %   >> radiusBall (Positive scalar)
            %      A positive scalar denoting the radius (in m) of the ball, R_0.
            %
            %   >> distancePivotPlaneBallContact (Positive scalar)
            %      A positive scalar denoting the distance (in m) from the pivot to the plane of ball contact on the beam, d.
            %
            %   >> distancePivotBeamCenter (Positive scalar)
            %      A positive scalar denoting the distance (in m) from the pivot to the center of mass of the beam, D.
            %      It holds D <= d.
            %
            %   >> distanceBallAxisBeam (Positive scalar)
            %      A positive scalar denoting the distance (in m) between the ball's center of gravity (its axis of rotiation)
            %      and its point of contact with the beam, R_1.
            %      It holds R_1 <= R_0.
            %
            %   >> frictionBeamBall (Nonnegative scalar)
            %      A positive scalar denoting the coefficient of viscous friction between the ball and the beam (in Ns/m).
            %
            %   >> frictionMotorBeam (Nonnegative scalar)
            %      A nonnegative scalar denoting the viscous damping/friction constant of the servo-motor (in Nm/(rad/s)).
            %
            %   >> electricalResistance (Positive scalar) 
            %      A positive scalar denoting the electrical resistance (in Ω) of the servo-motor.
            %
            %   >> electricalInductance (Positive scalar) 
            %      A positive scalar denoting the electrical inductance (in H) of the servo-motor.
            %
            %   >> motorTorqueConstant (Positive scalar) 
            %      A positive scalar denoting the electromotive torque constant (in Nm/A) of the servo-motor,
            %      equal to the back emf constant of the motor.
            %
            %   >> samplingInterval (Positive scalar, optional)
            %      The sampling interval (in seconds) that shall be used for the numerical simulation of the
            %      underlying ODE (the continuous-time plant model).
            %      If left out, the default vale 0.001 is used, which
            %      corresponds to a sampling rate of 1 kHz.         
            %
            % Returns:
            %   << this (BallBeamSystem)
            %      A new BallBeamSystem instance.
            
            this@NonlinearPlant(5, 1);
            
            this.massBeam = massBeam; % M
            this.lengthBeam = lengthBeam; % L
            this.inertiaBeam = inertiaBeam; % J_beam
            this.massBall = massBall; % m
            this.radiusBall = radiusBall; % R_0
            this.distancePivotPlaneBallContact = distancePivotPlaneBallContact; %d
            this.distancePivotBeamCenter = min(distancePivotBeamCenter, this.distancePivotPlaneBallContact); %D, D <= d
            this.distanceBallAxisBeam = min(distanceBallAxisBeam, this.radiusBall); %R_1, R_1 <= R_0
            this.frictionBeamBall = frictionBeamBall;
            this.frictionMotorBeam = frictionMotorBeam;
            this.electricalResistance = electricalResistance;
            this.electricalInductance = electricalInductance;
            this.motorTorqueConstant = motorTorqueConstant;
            
            if nargin > 13
                this.samplingInterval = samplingInterval;
            end         
        end
        
        %% setNoise
        function setNoise(this, noise)
            % Set the additional noise driving the plant, for simulation assumed constant over the sampling period.
            % The noise consists of two, typically
            % independent, components: A disturbance force (in newtons) driving the
            % movement of the ball, and an actuation noise (in volts) superimposing the control input (voltage) applied to the servo-motor.
            % Thus, the 2-dimensional joint distribution is passed here
            % with first component the marginal distribution of the
            % disturbance force acting on the ball. 
            %
            % Parameters:
            %   >> noise (Subclass of Distribution)
            %      The 2-dimensional system noise.
            
            assert(Checks.isClass(noise, 'Distribution') && noise.getDim() == 2, ...
                'DoubleInvertedPendulum:SetNoise:InvalidNoise', ...
                '** <noise> must be a 2-dimensional distribution **');                
            
            setNoise@SystemModel(this, noise);
                
            [~, noiseCov] = noise.getMeanAndCov();
            % store for convenience
            this.varDisturbanceForceActuator = noiseCov(2, 2);
            this.varDisturbanceForceBall = noiseCov(1, 1); 
        end
        
        %% linearizeAroundEquilibriumCont
        function [A_cont, B_cont, C_cont, G_cont, W_cont] = linearizeAroundEquilibriumCont(this, linContDisturbanceForceBall, linContDisturbanceForceAct)
            % Obtain a continuous-time linearization of the plant dynamics around the unstable equilibrium, where r=theta=i=0, 
            % that corresponds to the state [0 0 0 0 0]'.
            % The linearization is a differential equation of the form x_dot = A_cont * x + B_cont * u + G_cont * w 
            % and y = C_cont * x is the output.
            % The state x is given by [q, q_dot, i]' 
            % with q = [displacement of the ball's center of mass on the beam (r); angular rotation of the beam (theta)],
            % q_dot the time derivative of q, and i the current flowing through the servo-motor.
            % The input u is the voltage (in volts) applied to the servo-motor.
            % The noise w=[w1, w2]' is composed of a disturbance force (w1, in newtons) driving the movement of the ball 
            % and a disturbance (w2, in volts) acting on the servo-motor (additive on u, actuation noise). 
            % w1 and w2 are assumed white and mutually independent, the returned noise covariance W_cont is thus diagonal.
            %
            % Parameters:
            %   >> linContDisturbanceForceBall (Nonnegative scalar, optional)
            %      The variance of the disturbance force, in newtons
            %      squared, acting in the direction of the first
            %      generalized coordinate r, the ball position, thus
            %      driving the ball's movement.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForceBallContLin> property is used.
            %      If no disturbance torque is to be used, either set
            %      this property to 0 or pass 0 here.         
            %
            %   >> linContDisturbanceForceAct (Nonnegative scalar, optional)
            %      The variance of the disturbance, in volts squared,
            %      additionally entering at the actuator to drive the servo-motor.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForceActuatorContLin> property is used.
            %      If no disturbance force is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            % Returns:
            %   << A_cont (5-by-5 matrix)
            %      The system matrix of the continuous-time linearization.
            %
            %   << B_cont (5-dimensional column vector)
            %      The input matrix of the continuous-time linearization.
            %
            %   << C_cont (2-by-5 matrix)
            %      The measurement matrix, assuming that q as defined above
            %      is measured directly.
            %
            %   << G_cont (5-by-2 matrix)
            %      The system noise matrix, as described above, of the continuous-time linearization.
            %
            %   << W_cont (2-by-2 matrix)
            %      The covariance of the white noise acting on the
            %      continuous-time linearized dynamics as described above.
            %      Note: If either of the noise components is left out (variance = 0), this matrix is singular.
            
            arguments
                this
                linContDisturbanceForceBall (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceBallContLin;                
                linContDisturbanceForceAct (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceActuatorContLin;
            end
            
            mg =  BallBeamSystem.graviationalAcceleration * this.massBall;
            R1d = this.distanceBallAxisBeam + this.distancePivotPlaneBallContact;
            gMD = BallBeamSystem.graviationalAcceleration * this.massBeam * this.distancePivotBeamCenter;
            
            % we have M*q''+C*q'+G*q=tau+w for the linearized ball and beam system
            % without the electrical part
            % hence, as M is invertible, q''+M_inv*C*q'+M_inv*G*q=M_inv*tau+M_inv*w
            % mass matrix M, evaluated at r=0            
            M_lin11 = this.massBall * (1 + 0.4 * (this.radiusBall / this.distanceBallAxisBeam) ^ 2);
            M_lin12 = -this.massBall ...
                * (R1d + 0.4 * this.radiusBall ^ 2 / this.distanceBallAxisBeam);            
            M_lin22 = this.inertiaBeam + ...
                this.massBall * (R1d ^ 2 + 0.4 * this.radiusBall ^ 2); 
                         
            invM = [M_lin22, -M_lin12; -M_lin12, M_lin11] ./ (M_lin11 * M_lin22 - M_lin12 ^ 2);
            
            % compute based on linearized continuous-time dynamics
            A_cont = zeros(5);
            A_cont(1:2, 3:4) = eye(2);            
            A_cont(3:4, 1:2) = invM * [0, -mg; ...
                                        -mg, mg * R1d + gMD
                                        ];            
            A_cont(3:4, 3) = invM(:, 1) .* (-this.frictionBeamBall);
            A_cont(3:4, 4) = invM(:, 2) .* (-this.frictionMotorBeam);
            A_cont(3:4, 5) = invM(:, 2) .* this.motorTorqueConstant; 
            A_cont(5, 4) = -this.motorTorqueConstant / this.electricalInductance;
            A_cont(5, 5) = -this.electricalResistance / this.electricalInductance;            
            
            B_cont = [0; 0; 0; 0; 1/this.electricalInductance];            
      
            C_cont = eye(2, 5); % we measure r and theta
            
            % two noise components
            % first one acts an ball (r) directly
            % second one acts as disturbance input
            G_cont = zeros(5, 2);
            G_cont(3:4, 1) = invM(:, 1);
            G_cont(:, 2) = B_cont;

            W_cont = diag([linContDisturbanceForceBall; linContDisturbanceForceAct]); % semidefinite if one variance is 0
        end
        
         %% linearizeAroundEquilibrium
        function [A, B, C, W] = linearizeAroundEquilibrium(this, samplingInterval, contDisturbanceForceBall, contDisturbanceForceAct)
            % Obtain a discrete-time linearization of the dynamics around the unstable equilibrium 
            % that corresponds to the state [0 0 0 0 0]'.
            % The linearization is a difference equation of the form x_k = A * x_k-1 + B * u_k-1 + w_k-1 
            % and y_k = C * x_k is the output.
            % The state x_k is given by [q_k, q_k_dot, i_k]' 
            % with q_k = [displacement of the ball's center of mass on the beam (r_k); angular rotation of the beam (theta_k)]
            % and i_k the current flowing through the servo-motor.
            % The input u_k is the voltage (in volts) applied to the servo-motor, producing a torque to rotate the beam.
            % The noise w_k is the discrete-time equivalent of the continuous-time noise G_cont * w and has covariance
            % W = G_cont * W_cont * G_cont' with w, G_cont and W_cont as described and computed by linearizeAroundEquilibriumCont().
            % Note: If either of the noise components is left out (variance = 0), this matrix is singular.
            %
            % Parameters:
            %   >> samplingInterval (Positive scalar, optional)
            %      The sampling interval (in seconds) to be used for the
            %      discretization of the linearized continuous-time
            %      dynamics. If none is provided, the value given by this
            %      instance's <samplingInterval> property is used.
            %
            % Parameters:
            %   >> contDisturbanceForceBall (Nonnegative scalar, optional)
            %      The variance of the continuous-time disturbance force, in newtons
            %      squared, acting in the direction of the first
            %      generalized coordinate r, the ball position, thus
            %      driving the ball's movement.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForceBallContLin> property is used.
            %      If no disturbance torque is to be used, either set
            %      this property to 0 or pass 0 here.      
            %
            %   >> contDisturbanceForceAct (Nonnegative scalar, optional)
            %      The variance of the continuous-time disturbance, in volts squared,
            %      additionally entering at the actuator to drive the servo-motor.
            %      If none is provided, the value given by this
            %      instance's <varDisturbanceForceActuatorContLin> property is used.
            %      If no disturbance force is to be used, either set
            %      this property to 0 or pass 0 here.
            %
            % Returns:
            %   << A (5-by-5 matrix)
            %      The system matrix of the discrete-time linearization.
            %
            %   << B (5-dimensional column vector)
            %      The input matrix of the discrete-time linearization.
            %
            %   << C (2-by-5 matrix)
            %      The measurement matrix, assuming that q as defined above
            %      is measured directly at every sampling instant.
            %
            %   << W (2-by-2 matrix)
            %      The covariance matrix of the noise, as described above.
            
            arguments
                this                
                samplingInterval(1,1) double {mustBePositive} = this.samplingInterval;
                contDisturbanceForceBall (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceBallContLin;
                contDisturbanceForceAct (1,1) double {mustBeNonnegative, mustBeFinite} = this.varDisturbanceForceActuatorContLin;                
            end
            
            % compute based on linearized continuous-time dynamics
            [A_cont, B_cont, C_cont, G_cont, W_cont] = this.linearizeAroundEquilibriumCont(...
                contDisturbanceForceBall, contDisturbanceForceAct);
           
            % discretize A and B in one shot, assuming zero order hold for the
            % input (cf. Raymond DeCarlo, Linear Systems: A State Variable
            % Approach with Numerical Implementation, page 215)
            tmp = expm([A_cont B_cont; zeros(1,6)] * samplingInterval);
            A = tmp(1:5, 1:5);
            B = tmp(1:5, 6:end);
            C = C_cont;        
            % assume that noise is zero mean, white     
            G = G_cont * W_cont * G_cont';
            W = integral(@(x) expm(A_cont*x) * G * expm(A_cont'*x), 0, samplingInterval, 'ArrayValued', true);
            % symmetrize W
            W = (W + W') / 2;            
        end 
        
        %% isValidState
        function isValid = isValidState(this, state)
            % Function to check whether a given plant state is valid i.e.,
            % to check whether the ball has not reached the edn of the
            % beam.
            %
            % Parameters:
            %   >> state (5-dimensional vector)
            %      The system state to check.
            %
            % Returns:
            %   << isValid (Flag, i.e., boolean)
            %      Flag to indicate whether the given state state is valid.
            %      True is returned in case the ball position on the beam, r, is less than half the length of the beam, 
            %      and the beam's tilt angle, theta, is within [-45°,45°], false otherwise.
            %      More precisely, true is returned if |r| + R_0 < 0.5 * L and -45° <= theta <= 45°,
            %      where L is the length of the beam, R_0 is the radius of the ball,
            %      and theta is the angular rotation (tilt angle) of the beam.
            %      Otherwise false is returned.
            %      False is also returned if any of the state variables is +/-inf or +/-nan.
            
            arguments
                this
                state(5,1) double {mustBeReal} % if row vector is passed, it is converted automatically
            end            
            isValid = all(isfinite(state)) ...
                && 0.5 * this.lengthBeam > (abs(state(1)) + this.radiusBall) ...
                && 4 * abs(state(2)) <=  pi;
        end
    end
    
    methods (Access = protected)
        %% nonlinearDynamics
        function predictedStates = nonlinearDynamics(this, stateSamples, inputSamples, noiseSamples)
            tspan = [0 this.samplingInterval];                        
            mg = BallBeamSystem.graviationalAcceleration * this.massBall;
            MDg = BallBeamSystem.graviationalAcceleration * this.massBeam * this.distancePivotBeamCenter;
            R1d = this.distancePivotPlaneBallContact + this.distanceBallAxisBeam;
            mgR1d = mg * R1d;
                        
            % mass matrix is 2-by-2, only M(2,2) is state-dependent
            % M = [m11, m12; m12, m22]
            m11 = this.massBall * (1 + 0.4 * (this.radiusBall / this.distanceBallAxisBeam) ^ 2);
            m12 = -this.massBall * (R1d + 0.4 * this.radiusBall ^ 2 / this.distanceBallAxisBeam); %=m21
            detPortion = m12 ^ 2; % part of det(M)
            
            if isempty(inputSamples)
                inputSamples = zeros(1, size(stateSamples, 2));
            end
            
            if isempty(noiseSamples)
                noiseSamples = zeros(2, size(stateSamples, 2));
            end
            
            predictedStates = zeros(5, size(stateSamples, 2));
            for j=1:size(stateSamples, 2)                
                [~,x]=ode45(@(t, x) fun(t, x, inputSamples(:, j), noiseSamples(:, j)), tspan, stateSamples(:, j));
                predictedStates(:, j) = x(end, :)';
            end
            
            function dxdt = fun(~, x, u, w)
                % states in x are ordered as follows: q, q_dot, i, with q_dot = [r, theta]
                % and i (in amperes) the motor current
                
                % mass matrix M only depends on r=q(1)=x(1), last entry M(2,2)
                m22 = this.inertiaBeam + this.massBall * (R1d ^ 2 + x(1) ^ 2 + 0.4 * this.radiusBall^2);                
                % M = [m11, m12; m12, m22]
                % so it follows M^-1 = 1/det(M)*[m22,-m12; -m12, m11]
                
                % gravity terms, G(q)
                G = [mg * sin(x(2));
                     -sin(x(2)) * (mgR1d + MDg) + mg * x(1) * cos(x(2));
                 ];
                % Coriolis and centrifugal and damping terms,
                % C(q,q-dot)*q_dot
                C = [this.frictionBeamBall * x(3) - this.massBall * x(1) * (x(4) ^ 2);
                    2 * this.massBall * x(1) * x(3) * x(4) + this.frictionMotorBeam * x(4)];
                
                dxdt = zeros(5, 1);
                dxdt(1:2, :) = x(3:4); % q_dot=q_dot
                
                %d_dotdot =M(q)^-1*(-C(q,q_dot)*q_dot)-M(q)^-1*G(q)+M(q)^-1*[1;0]*w(1)+ M(q)^-1*[0;1]*tau
                % with applied torque tau=Km*i
                % M = [m11, m12; m12, m22]
                % so it follows M^-1 = 1/det(M)*[m22,-m12; -m12, m11]
                % and det(M) = m11*m22-m12^2
                dxdt(3:4) = ([m22, -m12; -m12, m11] * (-C - G + [w(1); this.motorTorqueConstant * x(5)])) ./ (m11 * m22 - detPortion);
                    
                % linear dynamics for i
                % u is input voltage
                dxdt(5) = (u + w(2) - this.electricalResistance * x(5) - this.motorTorqueConstant * x(4)) / this.electricalInductance;                
            end            
        end
    end    
    
end

