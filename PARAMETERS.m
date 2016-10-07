%------- SYSTEM= PARAMETERS

Bulk_beta=1053;  % 1/kT, Hartree^-1
density=0.0050;  % N/V,  Bohr^-3

%------- DISTANCE COARSING

Tape = inf;             % depricated
Level = 0.01;           % depricated

%------- SOLUTbE PARAMETERS

%energy_units = 'kcal';

%------- BULK PARAMETERS

Charge_o=-2*0.4238;     % Oxygen charge
Charge_h=0.4238;        % hydrogen charge

Sigma_o = 5.9828;                   % WaterOxygen diameter, Bohr
Epsilon_o = 2.476456137085393e-004; % WaterOxygen epsilon, Hartree

%Sigma_h = 1/0.5291772;
Sigma_h = 0.8 / 0.5291772;          % WaterHydrogen diameter, Bohr 
Epsilon_h = 7.3306e-005;            % WaterHydrogen epsilon, Hartree

r_oh=1.8897;         % Length of Water Oxygen-Hydrogen bond, Bohr
r_hh=3.0803;         % Distance between Water Hydrogens, Bohr


Write_Protocol_dGamma=0;


%------- CUTOFF PARAMETERS

Cutoff_t_ng=1.0;            % Cut-off parameter for Ng procedure
Cutoff_big_pot=1000;        % The cut-off top boundary of the potentials;
                            % We put the exp(-U) equals to zero if U>big_pot;
%------- BRIDGE PARAMETERS

LJ_SL_SV_Bridge_Mix = 0;        % Turn on/ Turn off bridge 

%------- CALC PARAMETERS

Lambda_Mix = 0.5;               % Lambda coupling parameter
%Lambda_CoarseSteps = 0.3;
Lambda_CoarseSteps = 0;
CoarseSteps_Treshold = 1;
Nested_Iterations =  0;

closure = @rHNC_closure;       % closure expression
                                % closure=@rMHNC_closure; % -> PL HNC
                                % closure=@rHNC_closure; % -> HNC

vdw_rules = @oplsaa_rules2;   % LJ mixing rules
                    % lorentz_barthelot_rules - > sigma+, eps*
                    % oplsaa_rules2 -> sigma*,eps*
% vdw_rulse: 
%   lorentz_barthelot_rules : sigma arithm, epsilon geom
%   oplsaa_rules2 : sigma geom, epsilon geom



potential_transformator = '';
MSACutFactor=1;
RCorr = 0;

Ng_function = @ng_fun_calc2;
Short_Range_FN = @short_range_pot2;
%-------
K_HK_File = 'K_HK.txt';
%SchemeGenerator = @NMGM_Scheme4x;
SchemeGenerator = 'NMGM_Scheme_25(1e-4);';
%SchemeGenerator = @()(NMGM_Scheme_25(1e-4));

NMGM_StartLevel = 0; 