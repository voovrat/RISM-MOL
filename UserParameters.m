%%% Units, used in this file:
%
%   Energy    - kcal/mol
%   Distance  - Angstroems
%   Charge    - positron charge
%   

%------- WATER PARAMETERS


% charges of water oxygen and hydrogen
user_HydrogenCharge = 0.4238;
user_OxygenCharge = -2*user_HydrogenCharge;

% diameters (LJ sigma) of water oxygen and hydrogen (Angstroems)
user_OxygenSigma = 3.166;
user_HydrogenSigma = 0.8;

% Lenneard Jones epsilon parameters of water oxygen and hydrogen ( kcal/mol)
user_OxygenEpsilon = 0.1554;
user_HydrogenEpsilon = 0.046;

% geometry of water molecule (distances between atoms, Angstroems)
user_DistanceOxygenHydrogen = 1;         % Length of Water Oxygen-Hydrogen bond, Angstr
user_DistanceHydrogenHydrogen = 1.63 ;   % Distance between Water Hydrogens, Angstr

%--------  MODEL

% **** Closure ****
% There are two possibilities for closure relation:
%
%       user_Closure = 'HNC'   - Hyper Netted Chain Closure    
%       user_Closure = 'PLHNC' - Partially Linearized Closure (default)
%
user_Closure='PLHNC';

% **** Lennard Jones Potential Mixing rules 
% Two different mixing rules are implemented to obtain pair potential
% parameters form atom parameters
% 
%  Geometric OPLSAA rules:
%               sigma_ab = sqrt(sigma_a*sigma_b)
%               epsilon_ab = sqrt(epsilon_a*epsilon_b)
%     to use these rules set: user_MixingRules = 'OPLSAA'; 
%
%  Lorentz-Berthelot rules :  
%               sigma_ab = (sigma_a+sigma_b)/2
%               epsilon_ab = sqrt( epsilon_a * epsilon_b )
%
%     to use these rules set:  user_MixingRules = 'LorentzBerthelot'; (default)
user_MixingRules = 'LorentzBerthelot';

% **** Repulsive Bridge Correction (indicates, whether or not Kovalenko
% repulsive bridge correction should be used.)
% To possible values: 'yes' or 'no'
% default: 'no'
user_UseRepulsiveBridge = 'no';

%------ NUMERICAL SCHEME

% Iteration Lambda Coupling Parameter (default 0.5)
user_LambdaCoupling = 0.5;

% Accuracy Threshold to stop the iterations ( default 1e-4 ) 
% (corresponds to 0.02 kcal/mol error in the free energy calculations
user_AccuracyThreshold = 1e-4;

