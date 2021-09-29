% load simulation networks
NetworkDelayDistribution = load('SimulatedNetworks.mat');
NetworkDelayDistribution = NetworkDelayDistribution.NetworkDelayDistribution;
NetworkDelayDistribution = NetworkDelayDistribution(2,:);

caDelayProb = NetworkDelayDistribution;
scDelayProb = NetworkDelayDistribution;
acDelayProb = NetworkDelayDistribution;
    
% System properties: discrete time double integrator plant (+noise) with default
% parameters 
% state: position of mass (1d), velocity of mass
% measurement: position is measured, but noisy
% control input: force applied to accelerate the mass
mass = 1; % kg
samplingInterval = 0.01; % sec, 100 Hz
A_cont = [0 1; 0 0];
B_cont = [0; 1 / mass];
W_cont = 0.1 * eye(2); % process noise, continuous time
C = [1 0]; % system is not observable if only velocity were observed

contSys = ss(A_cont, B_cont, C, 0); % no D matrix
discreteSys = c2d(contSys, samplingInterval);

W = integral(@(x) expm(A_cont*x) * W_cont * expm(A_cont'*x), ...
        0, samplingInterval, 'ArrayValued', true);

V = 0.2^2; % variance of the measurement noise

% for simplicity, use identiy matrices in the quadratic cost function
Q = eye(2); 
R = 1; 

initialPlantState =  [0; 0.2]; 
initialEstimate = Gaussian(initialPlantState, 0.5 * eye(2));

controlSequenceLength = 4;
maxMeasDelay = 6;

simTime = 1000; % number of timesteps to be simulated (10 seconds)

% draw the packet delays in advance
caPacketDelays = DiscreteSample(caDelayProb,simTime)-1;
scPacketDelays = DiscreteSample(scDelayProb,simTime)-1;
acPacketDelays = DiscreteSample(acDelayProb,simTime)-1;
    
ncs = NcsIntegratorExample(discreteSys.A, discreteSys.B, discreteSys.C, Q, R, W, V, ...
    simTime, caDelayProb, controlSequenceLength, maxMeasDelay);
results = ncs.simulate(initialPlantState, initialEstimate, caPacketDelays, scPacketDelays, acPacketDelays);

plot(results.realTrajectory(1,:))
%plot(zeros(1,simTime+1),'k')
