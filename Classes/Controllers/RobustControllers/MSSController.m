classdef MSSController < SequenceBasedController & ModelParamsChangeable & CaDelayProbsChangeable
    % Implementation of a mean square stabilizing linear controller for sequence-based control over
    % networks with time-varying and/or generally unknown delay and loss
    % probabilities.
    %
    % Literature:
    %   Florian Rosenthal and Uwe D. Hanebeck,
    %   Stability Analysis of Polytopic Markov Jump Linear Systems
    %   with Applications to Sequence-Based Control over Networks,
    %   IFAC PapersOnLine,
    %   Vol. 53, No. 2, pp. 3104-3111, 2020.
    
    % >> This function/class is part of CoCPN-Sim
    %
    %    For more information, see https://github.com/spp1914-cocpn/cocpn-sim
    %
    %    Copyright (C) 2019-2021  Florian Rosenthal <florian.rosenthal@kit.edu>
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
    
    properties (SetAccess = immutable, GetAccess = private)       
        F;
        G;
        dimEta;
    end
    
    properties (Access = private)
        % the input related part of the state, evolves according to 
        % eta_k+1=F*eta_k + G*U_k        
        etaState; % holds eta_{k-1}
        % matrices corresponding to the dynamics of the augmented system (polytopic MJLS)
        augA;
        augB;
        % the computed controller gain
        L;
        recomputeGain = false;
    end
       
    properties (SetAccess = private, GetAccess = ?MSSControllerTest)
        vertices; % vertices of the transition matrix polytope (called P in the paper)
        controllerDelta(1,1) double {mustBeNonnegative, mustBeLessThan(controllerDelta, 1)};
        
        % store already computed gains for (A,B) and delta
        % structure array
        % structure has fields A, B, L, and delta
        gains;
    end
    
    properties (SetAccess = immutable, GetAccess = public)
        assumeCorrelatedDelays (1,1) logical = false;
        useLmiLab(1,1) logical=false; 
        % by default, we do not use Matlab's internal solver LMI Lab from the
        % Robust Control Toolbox, but instead solve the feasibility problem
        % using yalmip and SDPT3
        % this is usually not only considerably faster and more reliable, but can produce slightly different results
    end
    
    properties (Access = private, Constant)
        SDPT3_MAXIT = 75; % default value is 50
        DELTA_CHAR = char(hex2dec('03b4')); % the greek letter delta (lower case)
        
        messages = containers.Map({
            'solverNotFound';
            'setup';
            'setupMss';
            'doneSetup'; 
            'callSolver'; 
            'doneCallSolver';
            'feasRepSolver';
            'infeasRepSolver';
            'numericalProbRepSolver';
            'lackProgRepSolver';
            'manCheckFeas';
            'doneManCheckFeas';
            'obtainGain';
            'doneObtainGain';
            'evalCond';
            'doneEvalCond';
            'computeBound';
            'doneComputeBound';
            'solverNotFound';
            'solverReportedError';
            'solverUnexpected';
            'solutionFeasible';
            'solutionInfeasible';
            'previouslyComputedGainAvail'
            }, ...
            {'** External solver (SDPT3) not found, falling back to internal solver (LMI Lab) **';
             '** Setting up feasibility problem **';
             '** Setting up feasibility problem to evaluate sufficient MSS condition **';
             '** Done: Setting up feasibility problem **';
             '** Calling external solver (SDPT3) **';
             '** Done: Calling external solver (SDPT3) **';
             '** Feasibility reported by external solver (SDPT3) **';
             '** Infeasibility reported by external solver (SDPT3) **';
             '** Numerical problems reported by external solver (SDPT3) **';
             '** Lack of progress reported by external solver (SDPT3) **';
             '** Checking feasibility of found solution manually **';
             '** Done: Checking feasibility of found solution manually **';
             '** Obtaining stabilizing gain (L) from %s solution **';
             '** Done: Obtaining stabilizing gain (L) from %s solution **';
             '** Lower bound is %f, evaluating sufficient MSS condition (Theorem 10) **';
             '** Done: Evaluating sufficient MSS condition (Theorem 10) **';
             '** Computing lower bound for joint spectral radius of set A_L **';
             '** Done: Computing lower bound for joint spectral radius of set A_L **';
             '** External solver (SDPT3) not found **';
             '** External solver (SDPT3) reported an error: %s **';
             '** External solver (SDPT3) reported an unexpected return code: %s **';
             '** Found solution seems to be feasible, consider calling isMeanSquareStable() **';
             '** Found solution seems to be infeasible, consider calling isMeanSquareStable() **';
             ['** Use previously computed gain (L) for (A,B) and ' MSSController.DELTA_CHAR '=%d **']
            });  
    end
    
    methods (Access = public)
        %% MSSController
        function this = MSSController(A, B, sequenceLength, delta, lazyInitGain, assumeCorrDelays, useLmiLab)
            % Class constructor.
            %
            % Parameters:
            %   >> A (Square Matrix)
            %      The system matrix of the plant.
            %
            %   >> B (Matrix)
            %      The input matrix of the plant.
            %            
            %   >> sequenceLength (Positive integer)
            %      N, the length of the input sequence (i.e., the number of
            %      control inputs) to be computed by the controller.
            %
            %   >> delta (Nonnegative scalar in [0,1))
            %      Upper bound for the last entry p_k,(N)(N) of all possible transition matrices P_k 
            %      of the augmented dynamical sytem (polytopic MJLS).
            %
            %   >> lazyInitGain (Flag, i.e., a logical scalar)
            %      Flag to indicate whether the computation of controller gain L shall be
            %      done during object instantiation, i.e., during this
            %      constructor call (lazyInitGain = false), or shall be deferred until the first control sequence is to be computed (lazyInitGain = true).
            %
            %   >> assumeCorrDelays (Flag, i.e., a logical scalar)
            %      Flag to indicate whether the vertices of the transition
            %      matrix polytope shall be computed under the assumption
            %      of Markovian delays and losses (with time-invariant
            %      transition probalities), or whether a white process
            %      shall be assumed (assumeCorrDelays = false).
            %      In the former case, the polytope has 2*N*(N+1)!
            %      vertices, in the latter case it has only 2*N
            %      vertices, with N the sequence length.
            %
            %   >> useLmiLab (Flag, i.e., a logical scalar, optional)
            %      Flag to indicate whether Matlab's LMI Lab shall be used 
            %      to solve the LMI feasibility problem required to
            %      synthesize a mean square stabilizing controller gain.
            %      Otherwise, yalmip and the infeasible path-following algorithm provided by 
            %      the open-source solver SDPT3 are used to set up and solve the LMIs which is typically
            %      considerably faster and more reliable.
            %      This, however, requires that SDPT3 is installed and
            %      accessible from Matlab.
            %      If left out, the default value false is used.
            %      If false is passed, and SDPT3 cannot be invoked, which
            %      is determined by checking if 'sdpt3' is available,
            %      this implementation false back to Matlab's LMI Lab.
            %      
            %      Literature:
            %       Reha H. Tütüncü, Kim-Chuan Toh, and Michael J. Todd,
            %       Solving semidefinite-quadratic-linear programs using SDPT3, 
            %       Mathematical Programming, Ser. B 95, pp. 189–217 (2003),
            %       https://doi.org/10.1007/s10107-002-0347-5.
            %       
            %       Johan Löfberg,
            %       YALMIP: a toolbox for modeling and optimization in MATLAB,
            %       Proceedings of the 2004 IEEE International Symposium on
            %       Computer Aided Control Systems Design,
            %       Taipei, Taiwan, 2004.
            %
            % Returns:
            %   << this (MSSController)
            %      A new MSSController instance.
            
            Validator.validateSystemMatrix(A);
            dimX = size(A,1);
            Validator.validateInputMatrix(B, dimX);
            dimU = size(B, 2);
            
            uncontrollableEigs = dimX - rank(ctrb(A, B));
            % Check for controllability of A by B
            assert(uncontrollableEigs == 0, ...
                'MSSController:InvalidPlant', ...
                '** Plant (A, B) has %d uncontrollable eigenvalue(s) **', uncontrollableEigs);
                            
            this = this@SequenceBasedController(dimX, dimU, sequenceLength, true);
            [this.F, this.G, ~, ~, this.augA, this.augB] = Utility.createAugmentedPlantModel(sequenceLength, A, B);
            this.dimEta = size(this.augA, 2) - dimX;
            % initially, the buffer is empty
            this.etaState = zeros(this.dimEta, 1);
            this.controllerDelta = delta;
            this.assumeCorrelatedDelays = assumeCorrDelays;
            
            if nargin > 6
                this.useLmiLab = useLmiLab;
            end
            if ~this.useLmiLab && exist('sdpt3', 'file') ~= 2
                % SDPT3 seems not to be on the path, fall back to internal
                % solver lmilab
                warning('MSSController:SolverNotFound', ...
                    MSSController.messages('solverNotFound'));
                this.useLmiLab = true;
            end
            
            this.initTransitionMatrixPolytope();
            if lazyInitGain
                this.recomputeGain = true;
            else
                this.computeAndSetGain();
            end
        end
        
        %% changeCaDelayProbs
        function changeCaDelayProbs(this, newCaDelayProbs)
            % first truncate to sequence length and
            % then check if last entry equals current controller delta
            % if not, recompute gain
            newProbs = Utility.truncateDiscreteProbabilityDistribution(newCaDelayProbs, this.sequenceLength + 1);
            if newProbs(end) ~= this.controllerDelta
                this.controllerDelta = newProbs(end);
                % the vertices of the polytope change (number of vertices remains)
                this.initTransitionMatrixPolytope();
                % remember to recompute gain L
                this.recomputeGain = true;
            end
        end
        
         %% changeModelParameters
        function changeModelParameters(this, newA, newB, ~)
            Validator.validateSystemMatrix(newA, this.dimPlantState);
            Validator.validateInputMatrix(newB, this.dimPlantState, this.dimPlantInput);
            
            [~, ~, H, ~] = Utility.createActuatorMatrices(this.sequenceLength, this.dimPlantInput);
            % in the the augmented model, F, G, H, J do not change
            % so only adapt the "first block row" in augA and augB            
            this.augA(1:this.dimPlantState, 1:this.dimPlantState, 1) = newA;
            for i=2:this.sequenceLength
                % update B*H for all but first and last mode
                % H is empty for first and last mode
                this.augA(1:this.dimPlantState, :, i) = [newA, newB * H(:, :, i)];
            end
            this.augA(1:this.dimPlantState, 1:this.dimPlantState, this.sequenceLength + 1) = newA;            
            % augB is only different for first mode (multiply B by J)
            % B*J is zero for all modes but first
            %
            this.augB(1:this.dimPlantState, 1:this.dimPlantInput, 1) = newB;
            
            this.recomputeGain = true;
        end
        
        %% getControllerGain
        function L = getControllerGain(this)
            % Get the parameter/gain of the linear, mode-independent control law.
            %
            % Returns:
            %
            %   << L (Matrix, might be empty)
            %      The controller gain L, such that U_k = L*state_k.
            %      The empty matrix is returned in case computation of gain 
            %      is deferred until the first control sequence U_0 is to be
            %      computed and this method is called prior to first
            %      invocation of computeControlSequence().
            %          
            L = this.L;
        end
        
        
        %% reset
        function reset(this)
            this.etaState = zeros(this.dimEta, 1);
        end
        
        %% setEtaState
        function setEtaState(this, newEta)
            assert(Checks.isVec(newEta, this.dimEta) && all(isfinite(newEta)), ...
                'MSSController:SetEtaState:InvalidEta', ...
                '** <newEta> must be a %d-dimensional vector **', ...
                this.dimEta);
            
            this.etaState = newEta(:);
        end
        
        %% isMeanSquareStable
        function [maxRho, isMSS] = isMeanSquareStable(this)
            % Try to determine whether the closed loop dynamics, 
            % given in terms of a polytopic Markov Jump Linear System,
            % is mean square stable. This is done based on Theorems 8 and
            % 10 given in the paper below.
            %    
            % Literature:
            %   Florian Rosenthal and Uwe D. Hanebeck,
            %   Stability Analysis of Polytopic Markov Jump Linear Systems 
            %   with Applications to Sequence-Based Control over Networks,
            %   IFAC PapersOnLine,
            %   Vol. 53, No. 2, pp. 3104-3111, 2020.
            %
            % Returns:
            %   << maxRho (Nonnegative scalar)
            %      The largest spectral radius of the matrices
            %      Lambda(0),..., Lambda(L) that form the set A_L associated
            %      with the second moment of the closed loop state,
            %      where L is the number of vertices of the polytope. This
            %      value is a natural lower bound of the joint spectral
            %      radius of A_L. If it is >= 1, then the closed loop
            %      dynamics is not mean square stable according to Theorem
            %      8 in the aforementioned paper, and hence <isMSS> is set to
            %      false.
            %
            %  << isMSS (Logical scalar, i.e., a flag)
            %     Flag to indicate whether the closed loop dynamics is mean
            %     square stable. 
            %     False is returned in case <maxRho> >= 1.
            %     If, however, false is returned and <maxRho> < 1, then the
            %     closed loop dynamics MAY OR MAY NOT be mean square
            %     stable, as the LMI condition provided by Theorem 10 in
            %     the aforementioned paper is only sufficient but not
            %     necessary. In such a case, <isMSS>=false should be considered a
            %     mere educated guess.
            %     True is returned if <maxRho> < 1 and the sufficient LMI
            %     condition by Theorem 10 in the paper is (numerically) feasible.
            
            % compute lower bound for the JSR of the set A_L associated with the
            % closed loop system: maximum of the spectral radii of each
            % matrix Lambda in A_L            
            disp(MSSController.messages('computeBound'));            
            numModes = this.sequenceLength + 1;
            dimState = size(this.augA, 1);                       
            A_cl = zeros(dimState, dimState, numModes);
            stackedA = [];
            for i=1:numModes
                A_cl(:, :, i) = this.augA(:, :, i) + this.augB(:, :, i) * this.L;
                stackedA = blkdiag(stackedA, kron(A_cl(:, :, i), A_cl(:, :, i)));
            end
            % compute the spectral radii
            I = eye(dimState^2);
            maxRho = 0;
            for r=1:size(this.vertices, 3)
                lambda = kron(this.vertices(:, :, r)', I) * stackedA;
                maxRho = max(max(abs(eig(lambda))), maxRho);
                if maxRho >= 1
                    break
                end
            end
            disp(MSSController.messages('doneComputeBound'));
            if maxRho >= 1
                % 1 <= maxRho <= JSR(A_L) -> closed loop dynamics not MSS
                % (Theorem 8)
                isMSS = false;                
                fprintf('** Lower bound is %f, by Theorem 8 the closed dynamics is not MSS **\n', maxRho);
            elseif this.useLmiLab                
                fprintf([MSSController.messages('evalCond') '\n'], maxRho);                
                isMSS = this.evaluateSufficientStabilityConditionLmiLab(A_cl);
                disp(MSSController.messages('doneEvalCond')); 
            else
                fprintf([MSSController.messages('evalCond') '\n'], maxRho);                
                isMSS = this.evaluateSufficientStabilityCondition(A_cl);
                disp(MSSController.messages('doneEvalCond')); 
            end
            % maxRho < 1 and isMSS = false does not necessarily imply that
            % closed loop dynamics is not MSS as LMI condition is only
            % sufficient and not necessary
            % maxRho < 1 and isMSS = true implies that
            % closed loop dynamics is MSS as LMI condition is sufficient
            % and was reported (numerically) feasible
            if maxRho < 1
                if isMSS
                     disp('** Closed loop dynamics is MSS **');
                else
                    disp('** Closed loop dynamics may or may not be MSS **');
                end
            end
        end
    end  
    
    methods (Access = protected)
        %% doControlSequenceComputation
        function inputSequence = doControlSequenceComputation(this, plantState, varargin)
            if this.recomputeGain
                this.computeAndSetGain();
                this.recomputeGain = false;
            end
            
            % varargin might be empty, we don't need it anyways
            % more general than in the paper: we get an estimate of the plant
            % state, which can of course be the true state
            [state, ~] = plantState.getMeanAndCov(); 
            
            inputSequence = this.L * [state; this.etaState];
            % update the eta portion of the augmented state
            this.etaState = this.F * this.etaState + this.G * inputSequence;
        end
        
        %% doStageCostsComputation
        function stageCosts = doStageCostsComputation(this, state, ~, ~)
            % input is not considered, simply compute norm of plant state
            stageCosts = norm(state);
        end
        
        %% doCostsComputation
        function costs = doCostsComputation(this, stateTrajectory, ~)
            % use vecnorm instead
            costs = arrayfun(@(k) norm(stateTrajectory(:, k)), 1:size(stateTrajectory, 2));
        end
    end
    
    methods (Access = private)        
        %% synthesizeController
        function res = synthesizeController(this)
            % Corollary 11 in the above-mentioned paper
            dimState = size(this.augA, 1);
            numModes = this.sequenceLength + 1;
            dimInput = size(this.augB, 2);
            numVertices = size(this.vertices, 3);
            I = eye(dimState);
            
            disp(MSSController.messages('setup'));
            % use yalmip in combination with SDPT3 (or Mosek, Sedumi)            
            E = sdpvar(dimState, dimState, 'symmetric'); % w.l.o.g. use symmetric E
            K = sdpvar(dimInput, dimState, 'full');
            Mir = cell(numVertices, numModes);
            Si = cell(1, numModes);
            AEBK = cell(1, numModes);            
            for i=1:numModes
                Si{i} = sdpvar(dimState, dimState, 'symmetric');
                AEBK{i} = (this.augA(:, :, i) * E + this.augB(:, :, i) * K)'; % (A(i)*E)^T + (B(i) * K)^T
            end       
            Sbar = blkdiag(Si{:});
            constraints = [];
            %constraints = [Sbar >= 0]; % implies that Si > 0
            for r=1:numVertices
                Pr = sqrt(squeeze(this.vertices(:, :, r))); % t(r) in the paper (18)
                for i=1:numModes                    
                    Z = AEBK{i} * kron(Pr(i, :), I);
                    Mir{r,i} = [2*E - Si{i}, Z; ... 
                        Z', Sbar];
                    constraints = [constraints, Mir{r,i} >= 0]; % implies Sbar > 0, and hence Si{i} > 0
                end
            end
            
            disp(MSSController.messages('doneSetup'));
            disp(MSSController.messages('callSolver')); 
            diagnostics = optimize(unblkdiag(constraints), [], MSSController.createYalmipOptsSdpt3());            
            disp(MSSController.messages('doneCallSolver'));            
            
            switch diagnostics.problem
                case 0
                    disp(MSSController.messages('feasRepSolver'));                                 
                    res.feasible = true;
                case {1, 4, 5} %1: deemed infeasible, 4: numerical problems, 5: lack of progress
                    if diagnostics.problem == 1
                         disp(MSSController.messages('infeasRepSolver')); 
                    elseif diagnostics.problem == 4
                        disp(MSSController.messages('numericalProbRepSolver')); 
                    else                        
                        disp(MSSController.messages('lackProgRepSolver')); 
                    end
                    disp(MSSController.messages('manCheckFeas'));                     
                    res.feasible = true;
                    % matrix E must be nonsingular, i.e., full rank
                    if rank(value(E)) ~= dimState
                        % matrix E seems to be singular (up to tolerance)
                        res.feasible = false;
                    else
                        % check if the symmetric matrices are all positive definite (up to tolerance)
                        for i = 1:numModes
                            for r=1:numVertices
                                eigVals = eig(value(Mir{r,i}));
                                isPosDef = all(eigVals > length(eigVals) * eps(max(eigVals)));
                                if ~isPosDef
                                    res.feasible = false;
                                    % infeasibility detected, so exit loop
                                    break;
                                end
                            end
                            if res.feasible
                                eigVals = eig(value(Si{i}));
                                isPosDef = all(eigVals > length(eigVals) * eps(max(eigVals)));
                                if ~isPosDef
                                    res.feasible = false;
                                    break
                                end
                            else
                                % infeasibility detected, so stop checking
                                break
                            end
                        end
                    end
                    % solution should be feasible
                    disp(MSSController.messages('doneManCheckFeas'));
                    if res.feasible
                        disp(MSSController.messages('solutionFeasible'));
                    else
                        disp(MSSController.messages('solutionInfeasible'));
                    end
                case -3
                    error('MSSController:SynthesizeController:SolverNotFound', ...
                        MSSController.messages('solverNotFound'));
                case 9
                    error('MSSController:SynthesizeController:ProblemInSolver', ...
                         MSSController.messages('solverReportedError'), diagnostics.info);
                otherwise      
                    error('MSSController:SynthesizeController:UnexpectedReturnCode', ...
                        MSSController.messages('solverUnexpected'), yalmiperror(diagnostics.problem));    
            end            
               
            if res.feasible
                msg{1} = sprintf(MSSController.messages('obtainGain'), 'feasible');
                msg{2} = sprintf(MSSController.messages('doneObtainGain'), 'feasible');
            else
                msg{1} = sprintf(MSSController.messages('obtainGain'), '(likely infeasible) found');
                msg{2} = sprintf(MSSController.messages('doneObtainGain'), '(likely infeasible) found');
            end
            disp(msg{1});
            res.gain = value(K) / value(E);            
            disp(msg{2});
        end
        
        %% synthesizeControllerLmiLab
        function res = synthesizeControllerLmiLab(this)
            % Corollary 11 in the above-mentioned paper
            dimState = size(this.augA, 1);
            numModes = this.sequenceLength + 1;
            dimInput = size(this.augB, 2);
            numVertices = size(this.vertices, 3);
            I = eye(dimState);

            disp(MSSController.messages('setup'));
            setlmis([]);
            E = lmivar(1, [dimState, 1]); % symmetric matrix w.o.l.g.
            K = lmivar(2, [dimInput, dimState]); % full matrix            
            % each Si is symmetric
            [Si, ~, sSi] = arrayfun(@(~) lmivar(1, [dimState 1]), 1:numModes, 'UniformOutput', false);
            Sbar = lmivar(3, blkdiag(sSi{:}));
            numLmis = 0;
            %lmiterm([-numLmis 1 1 Sbar], 1, 1); % 0 < Sbar, which implies 0 < Si
            for r=1:numVertices
                Pr = sqrt(squeeze(this.vertices(:, :, r)));
                % define LMIs of the form 0 < R
                for i=1:numModes              
                    numLmis = numLmis + 1;
                    lmiterm([-numLmis 1 1 E], 2, 1); % E+E' = 2E (E is symmetric)
                    lmiterm([-numLmis 1 1 Si{i}], -1, 1); % -Si{i}
                    lmiterm([-numLmis 2 2 Sbar], 1, 1);
                    lmiterm([-numLmis 1 2 -E], 1, this.augA(:, :, i)'* kron(Pr(i, :), I));
                    lmiterm([-numLmis 1 2 -K], 1, this.augB(:, :, i)'* kron(Pr(i, :), I));
                    % this lmiterm implies Sbar > 0 and thus Si{i} > 0
                 end
            end
            lmisys = getlmis;
            options = [0 0 0 0 1]; % mute the solver
            disp(MSSController.messages('doneSetup'));
            disp('** Calling LMI Lab (feasp) to solve LMI **');
            [tmin, xfeas] = feasp(lmisys, options); % tmin should be slightly negative for found solution to be strictly feasible
            disp('** Done: Calling LMI Lab (feasp) to solve LMI **');
            if tmin < 0
                disp('** Feasibility reported by LMI Lab (feasp) **');
                res.feasible = true;
            else
                disp('** Found solution not reported feasible by LMI Lab (feasp) **');
                disp(MSSController.messages('manCheckFeas'));
                % numerical problems reported, so check the feasibility manually
                res.feasible = true;                
                % matrix E must be nonsingular, i.e., full rank
                if rank(dec2mat(lmisys, xfeas, E)) ~= dimState
                    % matrix E seems to be singular (up to tolerance)
                    res.feasible = false;
                else
                    evals = evallmi(lmisys, xfeas);
                    for j=1:numLmis
                        [lhs, rhs] = showlmi(evals, j);
                        eigVals = eig(rhs-lhs); % must be positive definite
                        tol = length(eigVals) * eps(max(eigVals));
                        isPosDef = all(eigVals > tol);
                        if ~isPosDef % numerically sound check
                            res.feasible = false;                        
                            % infeasibility detected, so exit loop                          
                            break;
                        end
                    end
                end
                disp(MSSController.messages('doneManCheckFeas'));
                if res.feasible
                    % solution should be feasible                
                    disp(MSSController.messages('solutionFeasible'));
                else
                    disp(MSSController.messages('solutionInfeasible'));
                end
            end
            
            msg = {'** Obtaining stabilizing gain (L) from %s solution **';
                   '** Done: Obtaining stabilizing gain (L) from %s solution **'};     
            if res.feasible == true
                msg{1} = sprintf(msg{1}, 'feasible');
                msg{2} = sprintf(msg{2}, 'feasible');
            else
                msg{1} = sprintf(msg{1}, '(likely infeasible) found');
                msg{2} = sprintf(msg{2}, '(likely infeasible) found');
            end
            disp(msg{1});
            res.gain = dec2mat(lmisys, xfeas, K) / dec2mat(lmisys, xfeas, E);           
            disp(msg{2});
        end
        
        %% evaluateSufficientStabilityCondition
        function feasible = evaluateSufficientStabilityCondition(this, A_cl)
            % Theorem 10 in the paper, check if the closed loop system is
            % MSS (sufficient condition only)
            dimState = size(A_cl, 1);
            numModes = this.sequenceLength + 1;
            numVertices = size(this.vertices, 3);
            I = eye(dimState);
            
            disp(MSSController.messages('setupMss'));
            % use yalmip in combination with MOSEK (or Sedumi)
            Ei = cell(1, numModes);            
            Si = cell(1, numModes);            
            EA = cell(1, numModes);            
            Mir = cell(numVertices, numModes);
            for i=1:numModes                
                Ei{i} = sdpvar(dimState, dimState, 'symmetric'); % w.l.o.g. E can be symmetric
                Si{i} = sdpvar(dimState, dimState, 'symmetric');
                EA{i} = (A_cl(:, :, i) * Ei{i})'; % take the transpose of the product: (AE)^T
            end
            Sbar = blkdiag(Si{:});
            
            %constraints = [Sbar >= 0]; % implies that Si > 0
            constraints = [];
            for r=1:numVertices
                Pl = sqrt(squeeze(this.vertices(:, :, r))); % t(r) in the paper (18)
                for i=1:numModes                    
                    Z = EA{i} * kron(Pl(i, :), I);
                     Mir{r,i} = [2 * Ei{i} - Si{i}, Z; ... 
                        Z', Sbar];
                    constraints = [constraints, Mir{r,i} >= 0]; % Schur complement implies Sbar > 0, and thus also Si > 0
                end
            end
           
            disp(MSSController.messages('doneSetup'));            
            disp(MSSController.messages('callSolver'));                      
            diagnostics = optimize(unblkdiag(constraints), [], MSSController.createYalmipOptsSdpt3());
            disp(MSSController.messages('doneCallSolver'));
            switch diagnostics.problem
                case 0
                   disp(MSSController.messages('feasRepSolver'));               
                   feasible = true;
                case {1, 4, 5} %1: deemed infeasible, 4: numerical problems, 5: lack of progress
                    if diagnostics.problem == 1
                        disp(MSSController.messages('infeasRepSolver'));                        
                    elseif diagnostics.problem == 4
                        disp(MSSController.messages('numericalProbRepSolver'));                        
                    else
                        disp(MSSController.messages('lackProgRepSolver'));             
                    end
                    disp(MSSController.messages('manCheckFeas'));
                                        
                    feasMsg = '** Found solution seems to be feasible **';
                    infeasMsg = '** Found solution seems to be infeasible **';
                    
                    % numerical problems reported, so check the eigenvalues manually
                    feasible = true;
                    % matrices Ei must be nonsingular, i.e., have full rank
                    if any(cellfun(@(E) rank(value(E)), Ei) ~= dimState)
                        % at least one Ei seems to be singular (up to tolerance)
                        feasible = false;
                        disp(MSSController.messages('doneManCheckFeas'));                  
                        disp(infeasMsg);           
                        return
                    else
                        for i = 1:numModes                    
                            for r=1:numVertices
                                eigVals = eig(value(Mir{r,i}));
                                isPosDef = all(eigVals > length(eigVals) * eps(max(eigVals)));
                                if ~isPosDef
                                    feasible = false;
                                    disp(MSSController.messages('doneManCheckFeas'));
                                    disp(infeasMsg);
                                    return
                                end
                            end                            
                            eigVals = eig(value(Si{i}));
                            isPosDef = all(eigVals > length(eigVals) * eps(max(eigVals)));
                            if ~isPosDef
                                feasible = false;
                                disp(MSSController.messages('doneManCheckFeas'));
                                disp(infeasMsg);              
                                return
                            end
                        end
                    end
                    disp(MSSController.messages('doneManCheckFeas'));                   
                    disp(feasMsg);
                case -3
                    error('MSSController:EvaluateSufficientStabilityCondition:SolverNotFound', ...
                        MSSController.messages('solverNotFound'));
                case 9
                    error('MSSController:EvaluateSufficientStabilityCondition:ProblemInSolver', ...
                    MSSController.messages('solverReportedError'), diagnostics.info);
                otherwise
                    error('MSSController:EvaluateSufficientStabilityCondition:UnexpectedReturnCode', ...
                    MSSController.messages('solverUnexpected'), yalmiperror(diagnostics.problem));
            end            
        end
    
        %% evaluateSufficientStabilityConditionLmiLab
        function feasible = evaluateSufficientStabilityConditionLmiLab(this, A_cl)
            % Theorem 10 in the paper, check if the closed loop system is
            % MSS (sufficient condition only)
            dimState = size(A_cl, 1);
            numModes = this.sequenceLength + 1;
            numVertices = size(this.vertices, 3);
            I = eye(dimState);
            
            disp(MSSController.messages('setupMss'));
            setlmis([]);
            Si = cell(1, numModes);
            sSi = cell(1, numModes);
            Ei = cell(1, numModes);
            % Si and Ei are symmetric
            for i = 1:numModes
                [Si{i},~, sSi{i}] = lmivar(1, [dimState 1]);
                Ei{i} = lmivar(1, [dimState, 1]); % symmetric matrix w.o.l.g.
            end
            Sbar = lmivar(3, blkdiag(sSi{:}));
            numLmis = 0;
            %lmiterm([-numLmis 1 1 Sbar], 1, 1); % 0 < Sbar, which implies 0 < Si
            for r=1:numVertices
                Pr = sqrt(squeeze(this.vertices(:, :, r)));
                % define LMIs of the form 0 < R
                for i=1:numModes              
                    numLmis = numLmis + 1;
                    lmiterm([-numLmis 1 1 Ei{i}], 2, 1); % Ei{i}+Ei{i}' = 2Ei{i} (Ei{i} is symmetric)                    
                    lmiterm([-numLmis 1 1 Si{i}], -1, 1); % -Si{i}
                    lmiterm([-numLmis 2 2 Sbar], 1, 1);
                    lmiterm([-numLmis 1 2 -Ei{i}], 1, A_cl(:, :, i)'* kron(Pr(i, :), I));
                    % this lmiterm implies Sbar > 0 and thus Si{i} > 0
                 end
            end
            lmisys = getlmis;
            options = [0 0 0 0 1]; % mute the solver
            disp(MSSController.messages('doneSetup'));
            disp('** Calling LMI Lab (feasp) to evaluate sufficient MSS condition **');
            [tmin, xfeas] = feasp(lmisys, options); % tmin should be slightly negative for found solution to be strictly feasible
            disp('** Done: Calling feasp to evaluate sufficient MSS condition **');
            if tmin < 0
                disp('** Feasibility reported by LMI Lab (feasp) **');
                feasible = true;
            else
                disp('** Found solution not reported feasible by LMI Lab (feasp) **');
                disp(MSSController.messages('manCheckFeas'));
                             
                % numerical problems reported, so check the feasibility manually
                feasible = true;
                % matrices Ei must be nonsingular, i.e., full rank
                ranks = cellfun(@(E) rank(dec2mat(lmisys, xfeas, E)), Ei);
                if any(ranks) ~= dimState
                    % at least one Ei seems to be singular (up to tolerance)
                    feasible = false;
                else
                    evals = evallmi(lmisys, xfeas);
                    for j=1:numLmis
                        [lhs, rhs] = showlmi(evals, j);
                        eigVals = eig(rhs-lhs); % must be positive definite
                        tol = length(eigVals) * eps(max(eigVals));
                        isPosDef = all(eigVals > tol);
                        if ~isPosDef % numerically sound check
                            feasible = false;                      
                            % infeasibility detected, so exit loop                          
                            break;
                        end
                    end
                end
                disp(MSSController.messages('doneManCheckFeas'));   
                if feasible
                    disp('** Found solution seems to be feasible **');
                else
                    disp('** Found solution seems to be infeasible **');
                end
            end           
        end
        
        
        
        %% initTransitionMatrixPolytope
        function initTransitionMatrixPolytope(this)
            if this.assumeCorrelatedDelays
                % we use the multi-simplex approach since we cannot exploit
                % any dependencies between the transition probailities
                % we only know that possible matrices are lower Hessenberg
                % matrices and that last entry is bounded from above by
                % delta
                % hence, general structure of P_k is:
                %
                % [? ? 0 0 0
                %  ? ? ? 0 0
                %  ? ? ? ? 0
                %  ? ? ? ? ?
                %  ? ? ? ? [0,delta]]
                %
                % where ? means that the corresponding entry is arbitrary
                % and only known to lie in [0,1]

                numModes = this.sequenceLength + 1;                
                vertsPerRow = cell(1, numModes);
                idxPerRow = cell(1, numModes);
                for j=1:this.sequenceLength
                    % j+1 vertices for the i-th row
                    % first row: 1-simplex, so 2 vertices
                    % second row: 2-simplex, so 3 vertices, and so forth
                    vertsPerRow{j} = eye(j+1, numModes);
                    idxPerRow{j} = 1:j+1;
                end
                if this.controllerDelta == 0
                    numVertices = this.sequenceLength * factorial(numModes);
                    % for the last row, we have sequenceLength many vertices with 1 at position i
                    vertsPerRow{numModes} = eye(this.sequenceLength, numModes);
                    idxPerRow{numModes} = 1:this.sequenceLength;
                else
                    numVertices = 2 * this.sequenceLength * factorial(numModes); % 2*N*(N+1)! vertices
                    % for the last row, we have 2*sequenceLength vertices
                    % sequenceLength many with 1 at position i
                    % sequenceLength many with 1-delta at position i and delta
                    % at last position
                    vertsPerRow{numModes} = [  eye(this.sequenceLength, numModes);
                                               [(1-this.controllerDelta) * eye(this.sequenceLength), repmat(this.controllerDelta, this.sequenceLength, 1)]
                                    ];
                    idxPerRow{numModes} = 1:2*this.sequenceLength;
                end
               
                [T{1:numModes}] = ndgrid(idxPerRow{:});
                idx = cell2mat(cellfun(@(Ti) Ti(:), T, 'UniformOutput', false));                
                % simply enumerate all vertices by combining their
                % respective indices
                this.vertices = zeros(numModes, numModes, numVertices);
                for j = 1:numVertices
                    for i=1:numModes
                        this.vertices(i, :, j) = vertsPerRow{i}(idx(j, i), :);
                    end
                end
            else
                this.computePolytopeWhiteDelays();                
            end
        end

        %% computePolytopeWhiteDelays
        function computePolytopeWhiteDelays(this)
           % for the case of white delays: exploit the Hessenberg
            % structure and the element interdependencies
            % assume delta is [0,1]
                
            switch this.sequenceLength
                case 1 % sequence length is 1
                    % 2 vertices
                    this.vertices(:, :, 2) = [1-this.controllerDelta, this.controllerDelta;
                                              1-this.controllerDelta, this.controllerDelta];
                    this.vertices(:,:, 1) = [1, 0;
                                             1, 0];
                case 2 % sequence length is 2
                    % 4=2+1+1 vertices
                    % including 2 for corner case of delay probs Pr[delay = 0] = 1 and Pr[delay = 1] = 1
                    this.vertices(:, :, 4) = [1 0 0;
                                              1 0 0;
                                              1 0 0];
                    this.vertices(:, :, 3) = [0 1 0;
                                              0 1 0;
                                              0 1 0];                    
                    this.vertices(:, :, 2) = [1-this.controllerDelta, this.controllerDelta, 0;
                                              1-this.controllerDelta, 0, this.controllerDelta;
                                              1-this.controllerDelta, 0, this.controllerDelta];
                    this.vertices(:, :, 1) = [0, 1, 0;
                                              0, 1-this.controllerDelta, this.controllerDelta;
                                              0, 1-this.controllerDelta, this.controllerDelta];
                case 3 % sequence length is 3, 6 vertices
                    this.initTransitionMatrixPolytopeSeq3(); 
                case 4 % sequence length is 4, 8 vertices
                    this.initTransitionMatrixPolytopeSeq4();
                case 5 % sequence length is 5, 10 vertices
                    this.initTransitionMatrixPolytopeSeq5();
                otherwise
                    % general case: we have 2*N vertices for the last row
                    % sequenceLength many with 1 at position i
                    % sequenceLength many with 1-delta at position i and delta
                    % at last position
                    vertsLastRow = [eye(this.sequenceLength, this.sequenceLength + 1);
                                 [(1-this.controllerDelta) * eye(this.sequenceLength), repmat(this.controllerDelta, this.sequenceLength, 1)]
                                ];
                    for j=1:2*this.sequenceLength
                        sums = cumsum(vertsLastRow(j, :));
                        this.vertices(:, :, j) = diag(1 - sums(1:end -1), 1);
                        % adapted from Utility.calculateDelayTransitionMatrix
                        % to create transition matrix without normalization
                        % of probabilities
                        for i=1:this.sequenceLength + 1
                            this.vertices(i:end, i, j) = vertsLastRow(j, i);
                        end
                    end                    
            end
            if this.controllerDelta == 0
                % remove redundant slices, we have way fewer slices 
                % (only N, corresponding to the corner cases mentioned above)                
                tmp = reshape(this.vertices, this.sequenceLength + 1, [], 1);
                tmp = reshape(tmp(:), (this.sequenceLength + 1) ^ 2, []);
                tmp = unique(tmp', 'rows', 'stable');
                this.vertices = reshape(tmp', this.sequenceLength + 1, this.sequenceLength + 1, []);
                if this.sequenceLength > 5
                    % get rid of vertex with 1 at last position (if present, then exactly one vertex)
                    for j=1:size(this.vertices, 3)
                        this.vertices(end, end, j)
                        if this.vertices(end, end, j) == 1
                            this.vertices = cat(3, this.vertices(:, :, 1:j-1), this.vertices(:, :, j+1:end));
                            break
                        end
                    end
                end
            end
        end

        %% initTransitionMatrixPolytopeSeq3
        function initTransitionMatrixPolytopeSeq3(this)
            % 6 vertices
            % including 3 for corner cases of delay probs Pr[delay = 0] = 1
            % Pr[delay = 1] = 1 and Pr[delay = 2] = 1            
            
            this.vertices(:, :, 6) = [1 0 0 0;
                                      1 0 0 0;
                                      1 0 0 0;
                                      1 0 0 0];
            
            this.vertices(:, :, 5) = [0 1 0 0;
                                      0 1 0 0;
                                      0 1 0 0;
                                      0 1 0 0];
            
            this.vertices(:, :, 4) = [0 1 0 0;
                                     0 0 1 0;
                                     0 0 1 0;
                                     0 0 1 0];
            %
            this.vertices(:, :, 3) = [1-this.controllerDelta, this.controllerDelta, 0, 0;
                                      1-this.controllerDelta, 0, this.controllerDelta, 0;
                                      1-this.controllerDelta, 0, 0, this.controllerDelta;
                                      1-this.controllerDelta, 0, 0, this.controllerDelta];
            this.vertices(:, :, 2) = [0, 1, 0, 0;
                                      0, 1-this.controllerDelta, this.controllerDelta, 0;
                                      0, 1-this.controllerDelta, 0, this.controllerDelta;
                                      0, 1-this.controllerDelta, 0, this.controllerDelta];
            this.vertices(:, :, 1) = [0, 1, 0, 0;
                                      0, 0, 1, 0;
                                      0, 0, 1-this.controllerDelta, this.controllerDelta;
                                      0, 0, 1-this.controllerDelta, this.controllerDelta];
            %
        end
        
        %% initTransitionMatrixPolytopeSeq4
        function initTransitionMatrixPolytopeSeq4(this)
            % 8
            % incl. 4 for corner cases of delay probs Pr[delay = 0] = 1
            % Pr[delay = 1] = 1, Pr[delay = 2] = 1 and Pr[delay = 3] = 1
            % so 2*N in total, where N=4 (sequence length)
            
            this.vertices(:, :, 8) = [1 0 0 0 0;
                                      1 0 0 0 0;
                                      1 0 0 0 0;
                                      1 0 0 0 0;
                                      1 0 0 0 0];           
            this.vertices(:, :, 7) = [0 1 0 0 0;
                                      0 1 0 0 0;
                                      0 1 0 0 0;
                                      0 1 0 0 0;
                                      0 1 0 0 0];
            %
            this.vertices(:, :, 6) = [0 1 0 0 0;
                                      0 0 1 0 0;
                                      0 0 1 0 0;
                                      0 0 1 0 0;
                                      0 0 1 0 0];
            this.vertices(:, :, 5) = [0 1 0 0 0;
                                      0 0 1 0 0;
                                      0 0 0 1 0;
                                      0 0 0 1 0;
                                      0 0 0 1 0];
            %
            this.vertices(:, :, 4) = [1-this.controllerDelta, this.controllerDelta, 0, 0, 0;
                                      1-this.controllerDelta, 0, this.controllerDelta, 0, 0;
                                      1-this.controllerDelta, 0, 0, this.controllerDelta, 0;
                                      1-this.controllerDelta, 0,0,0, this.controllerDelta;
                                      1-this.controllerDelta, 0,0,0, this.controllerDelta];
            this.vertices(:, :, 3) = [0, 1, 0, 0, 0;
                                      0, 1-this.controllerDelta, this.controllerDelta, 0, 0;
                                      0, 1-this.controllerDelta, 0, this.controllerDelta, 0;
                                      0, 1-this.controllerDelta, 0, 0, this.controllerDelta;
                                      0, 1-this.controllerDelta, 0, 0, this.controllerDelta];
            %
            this.vertices(:, :, 2) = [0, 1, 0, 0, 0;
                                      0, 0, 1, 0, 0;
                                      0, 0, 1-this.controllerDelta, this.controllerDelta, 0;
                                      0, 0, 1-this.controllerDelta, 0, this.controllerDelta;
                                      0, 0, 1-this.controllerDelta, 0, this.controllerDelta];
            this.vertices(:, :, 1) = [0, 1, 0, 0, 0;
                                      0, 0, 1, 0, 0;
                                      0, 0, 0, 1, 0;
                                      0, 0, 0, 1-this.controllerDelta, this.controllerDelta;
                                      0, 0, 0, 1-this.controllerDelta, this.controllerDelta];
        end
        
        %% initTransitionMatrixPolytopeSeq5
        function initTransitionMatrixPolytopeSeq5(this)
            % 10 vertices
            % incl. 5 for corner cases of delay probs Pr[delay = 0] = 1
            % Pr[delay = 1] = 1, Pr[delay = 2] = 1, Pr[delay = 3] = 1 and Pr[delay = 4] = 1
            % so 2*N in total, where N=5 (sequence length)
            this.vertices(:, :, 10) = [1 0 0 0 0 0;
                                       1 0 0 0 0 0;
                                       1 0 0 0 0 0;
                                       1 0 0 0 0 0;
                                       1 0 0 0 0 0;
                                       1 0 0 0 0 0];           
            this.vertices(:, :, 9) = [0 1 0 0 0 0;
                                       0 1 0 0 0 0;
                                       0 1 0 0 0 0;
                                       0 1 0 0 0 0;
                                       0 1 0 0 0 0;
                                       0 1 0 0 0 0];
            %
            this.vertices(:, :, 8) = [0 1 0 0 0 0;
                                       0 0 1 0 0 0;
                                       0 0 1 0 0 0;
                                       0 0 1 0 0 0;
                                       0 0 1 0 0 0;
                                       0 0 1 0 0 0];
            this.vertices(:, :, 7) = [0 1 0 0 0 0;
                                       0 0 1 0 0 0;
                                       0 0 0 1 0 0;
                                       0 0 0 1 0 0;
                                       0 0 0 1 0 0;
                                       0 0 0 1 0 0];
            %
            this.vertices(:, :, 6) = [0 1 0 0 0 0;
                                       0 0 1 0 0 0;
                                       0 0 0 1 0 0;
                                       0 0 0 0 1 0;
                                       0 0 0 0 1 0;
                                       0 0 0 0 1 0];      
            this.vertices(:, :, 5) = [  1-this.controllerDelta, this.controllerDelta, 0, 0, 0, 0;
                                        1-this.controllerDelta, 0, this.controllerDelta, 0, 0, 0;
                                        1-this.controllerDelta, 0, 0, this.controllerDelta, 0, 0;
                                        1-this.controllerDelta, 0, 0, 0, this.controllerDelta, 0;
                                        1-this.controllerDelta, 0, 0, 0, 0, this.controllerDelta;
                                        1-this.controllerDelta, 0, 0, 0, 0, this.controllerDelta];
            %
            this.vertices(:, :, 4) = [  0, 1, 0, 0, 0, 0;
                                        0, 1-this.controllerDelta, this.controllerDelta, 0, 0, 0;
                                        0, 1-this.controllerDelta, 0, this.controllerDelta, 0, 0;
                                        0, 1-this.controllerDelta, 0, 0, this.controllerDelta, 0;
                                        0, 1-this.controllerDelta, 0, 0, 0, this.controllerDelta;
                                        0, 1-this.controllerDelta, 0, 0, 0, this.controllerDelta];
            this.vertices(:, :, 3) = [  0, 1, 0, 0, 0, 0;
                                        0, 0, 1, 0, 0, 0;
                                        0, 0, 1-this.controllerDelta, this.controllerDelta, 0, 0;
                                        0, 0, 1-this.controllerDelta, 0, this.controllerDelta, 0;
                                        0, 0, 1-this.controllerDelta, 0, 0, this.controllerDelta;
                                        0, 0, 1-this.controllerDelta, 0, 0, this.controllerDelta];
            %
            this.vertices(:, :, 2) = [  0, 1, 0, 0, 0, 0;
                                        0, 0, 1, 0, 0, 0;
                                        0, 0, 0, 1, 0, 0;
                                        0, 0, 0, 1-this.controllerDelta, this.controllerDelta, 0;
                                        0, 0, 0, 1-this.controllerDelta, 0, this.controllerDelta;
                                        0, 0, 0, 1-this.controllerDelta, 0, this.controllerDelta];
            this.vertices(:, :, 1) = [  0, 1, 0, 0, 0, 0;
                                        0, 0, 1, 0, 0, 0;
                                        0, 0, 0, 1, 0, 0;
                                        0, 0, 0, 0, 1, 0;
                                        0, 0, 0, 0, 1-this.controllerDelta, this.controllerDelta;
                                        0, 0, 0, 0, 1-this.controllerDelta, this.controllerDelta];
        end
        
        %% computeAndSetGain
        function computeAndSetGain(this)
            currA = this.augA(1:this.dimPlantState, 1:this.dimPlantState, this.sequenceLength + 1); % A
            currB = this.augB(1:this.dimPlantState, 1:this.dimPlantInput, 1); % B
            % search, if gains had been computed previously
            for i=1:numel(this.gains)
                if this.gains(i).delta == this.controllerDelta && isequal(this.gains(i).A, currA) && isequal(this.gains(i).B, currB)
                    this.L = this.gains(i).L;
                    fprintf([MSSController.messages('previouslyComputedGainAvail') '\n'], this.controllerDelta);
                    return
                end
            end

            if this.useLmiLab
                res = this.synthesizeControllerLmiLab();
            else
                res = this.synthesizeController();
            end
            this.L = res.gain;

            % store information, i.e., append to struct array
            this.gains(end + 1).A = currA;
            this.gains(end).B = currB;
            this.gains(end).L = this.L;
            this.gains(end).delta = this.controllerDelta;
        end
    end
    
    methods (Access = private, Static)
        %% createYalmipOptsSdpt3
        function opts = createYalmipOptsSdpt3()
            opts = sdpsettings('solver', 'sdpt3');
            
            opts.verbose = 0; % mute the solver
            opts.sdpt3.maxit = MSSController.SDPT3_MAXIT;
        end
    end
end

