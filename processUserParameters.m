%%% Units, used in this file:
%
%   Energy    - kcal/mol
%   Distance  - Angstroems
%   Charge    - positron charge
%   

[DU,EU]=unit_names;

%------- WATER PARAMETERS
% charges of water oxygen and hydrogen
if(ismember('user_OxygenCharge',who))
    Charge_o = user_OxygenCharge;
end
   
if(ismember('user_HydrogenCharge',who)) 
    Charge_h = user_HydrogenCharge;
end


% diameters (LJ sigma) of water oxygen and hydrogen (Angstroems)

if(ismember('user_OxygenSigma',who))
    Sigma_o = user_OxygenSigma * unit2unit(DU,'Angstr','Bohr');
end
if(ismember('user_HydrogenSigma',who))
    Sigma_h = user_HydrogenSigma * unit2unit(DU,'Angstr','Bohr');
end

% Lenneard Jones epsilon parameters of water oxygen and hydrogen ( kcal/mol)

if(ismember('user_OxygenEpsilon',who))
    Epsilon_o = user_OxygenEpsilon * unit2unit(EU,'kcal/mol','Hartree');
end
if(ismember('user_HydrogenEpsilon',who))
    Epsilon_h = user_HydrogenEpsilon * unit2unit(EU,'kcal/mol','Hartree');
end

% geometry of water molecule (distances between atoms, Angstroems)
if(ismember('user_DistanceOxygenHydrogen',who))
    r_oh = user_DistanceOxygenHydrogen *  unit2unit(DU,'Angstr','Bohr');
end
if(ismember('user_DistanceHydrogenHydrogen',who))
    r_hh = user_DistanceHydrogenHydrogen *  unit2unit(DU,'Angstr','Bohr');
end

%--------  MODEL

% **** Closure ****
% There are two possibilities for closure relation:
%
%       user_Closure = 'HNC'   - Hyper Netted Chain Closure    
%       user_Closure = 'PLHNC' - Partially Linearized Closure (default)
%
if ismember('user_Closure',who)
    if strcmp(user_Closure,'HNC')
        closure = @rHNC_closure;
    elseif strcmp(user_Closure,'PLHNC')
        closure = @rMHNC_closure;
    else
        error([ 'Unknown closure: ' user_Closure ' (should be either ''HNC'' or ''PLHNC'' )']);
    end
end

% **** Lennard Jones Potential Mixing rules 
% Two different mixing rules are implemented to obtain pair potential
% parameters form atom parameters
% 
%  Geometric OPLSAA rules:
%               sigma_ab = sqrt(sigma_a*sigma_b)
%               epsilon_ab = sqrt(epsilon_a*epsilon_b)
%
%  Lorentz-Berthelot rules :  
%               sigma_ab = (sigma_a+sigma_b)/2
%               epsilon_ab = sqrt( epsilon_a * epsilon_b )
%
%     to use these rules:  user_MixingRules = 'LorentzBerthelot'; (default)


if ismember('user_MixingRules',who)
       
    if strcmp(user_MixingRules,'LorentzBerthelot')
        vdw_rules = @lorentz_barthelot_rules;
    elseif strcmp(user_MixingRules,'OPLSAA')
        vdw_rules = @oplsaa_rules2;
    else
        error([ 'Unknown mixing rules: ' uuser_MixingRules ' (should be either ''LorentzBerthelot'' or ''OPLSAA'' )']);
    end
end

% **** Repulsive Bridge Correction (indicates, whether or not Kovalenko
% repulsive bridge correction should be used.)
% To possible values: 'yes' or 'no'
% default: 'no'
if ismember('user_UseRepulsiveBridge',who)
    RBYN = upper(user_UseRepulsiveBridge);
    
    if strcmp(RBYN,'YES')
        LJ_SL_SV_Bridge_Mix = 1;
    elseif strcmp(RBYN,'NO')
        LJ_SL_SV_Bridge_Mix = 0;
    else
        error([ 'Unknown repulsive bridge switch: ' user_UseRepulsiveBridge ' (should be either ''yes'' or ''no'' )']);
    end
end

%------ NUMERICAL SCHEME

% Iteration Lambda Coupling Parameter (default 0.5)
if ismember('user_LambdaCoupling',who)
    Lambda_Mix = user_LambdaCoupling;
end

% Accuracy Threshold to stop the iterations ( default 1e-4 ) 
% (corresponds to 0.02 kcal/mol error in the free energy calculations
if ismember('user_AccuracyThreshold',who)
    SchemeGenerator = ['NMGM_Scheme_25(' num2str(user_AccuracyThreshold) ');'];
%    SchemeGenerator = @()(NMGM_Scheme_25(user_AccuracyThreshold));
end
