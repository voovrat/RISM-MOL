function [Mu,g,grid,gamma,c,pot_long,Err,timing,LJ_SL_SV,EXPONENTA,u_short,F_S,F_L,F_L_K] = startRISM(fname,input_units,output_units,xparam)
%
% usage: [Mu,g,grid,gamma,c] = startRISM(fname,input_units,output_units,xparam)
%
% INPUT PARAMETERS:
%
%   fname -     path to file, contining structure and charges of solute
%               molecule
%               Format of file: six columns of numbers (wiuthout a header),
%               separated by spaces or tabs
%
%               Columns are: X Y Z sigma epsilon charge
%                    X Y Z           - coordinates of atom 
%                    sigma   - sigma parameter of Van-der-Waals potential
%                              
%                    epsilon - epsilon parameter of Van-der-Waals potential
%                               
%                    charge   - charge of atom (atomic units)
%
%    input_units - structure, determines in which units is input file.
%    output_units - structure, determines in which units will be output
%    parameters ( Mu and grid )
%    
%        units.Energy - energy units.
%               supported: 'Hartree', 'J','KBT' (298K) ,'eV','kJ/mol','kcal/mol'
%
%        units.Distance - distance units.
%               supported: 'Angstr', 'Bohr', 'm'
%
%        by default: units.Energy = 'kcal/mol', units.Distance ='Angstr'
%
%
%   xparam - extra parameters. Only for experienced users.
%
% OUTPUT PARAMETERS:
%
%   Mu - styructure, containing chemical potentials calculated by different
%   methods
%   Fields of structure are: 
%       Mu.HNC - Hypper-Netted-Chain  (kcal/mol)
%       Mu.G   - Gaussian Fluctuations (kcal/mtypeol)
%       Mu.HNCB - Hypper-Netted-Chain with Repulsive Bridge correction (kcal/mol)
%       Mu.PW   - Partial Wave  (kcal/mol)
%       Mu.PWC  - Partial Wave Corrected (kcal/mol)
%       Mu.VUA  - Excluded Volume (Bohr^-3)
%
%   g - RDF functions array. 
%       The odd columns correspons to interaction between
%       solute atoms and water Oxygen
%       The even columns corresponds to interaction between solute atoms
%       and the water Hydrogen
%
%   grid -  grid vectors for RDFs 
%
%   gamma, c - solutions of the RISM equations. (short part)
%       RDF functions g = gamma + c + 1
%
%   pot_long - long-term potential of electristatic interactions
%


%[DI,EI,DTab,ETab] = build_unit_index;


[DU,EU]=unit_names;


if nargin<2
    input_units = struct('Energy','kcal/mol','Distance','Angstr');
end

if nargin<3
    output_units = struct('Energy','kcal/mol','Distance','Angstr');
end

if nargin<4
    xparam='';
end

xparam = [ 'UserParameters;' xparam '; processUserParameters;' xparam ];

silent = 0;
eval(xparam);

if ~silent
    fprintf('Input Units:\n');
    disp(input_units);

    fprintf('Output Units:\n');
    disp(output_units);
end
   
%K_HK_File = 'K_HK.txt';
eval(['PARAMETERS;' xparam ]);

if ismember('Initial_Guess_File',who)
    rgamma0 =  load(Initial_Guess_File);
    R0 = rgamma0(:,1);
    gamma0 = rgamma0(:,2:end);
else
    R0=NaN;
    gamma0=NaN;

end

t0=mytic;
[gamma,c,Err,timing,LJ_SL_SV,EXPONENTA,u_short,F_S,F_L,F_L_K]=run_NMGM_iterations(fname,input_units,K_HK_File,['PARAMETERS;' xparam],eval(SchemeGenerator),R0,gamma0);

if ~silent
    dt =mytoc(t0)
end
  timing.t_total=mytoc(t0);


if sum(sum(isnan(gamma)))
    Mu=NaN;
    g=NaN;
    grid=NaN;
    c=NaN;
    pot_long=NaN;
    return
end

g= gamma + c +1;

%gridBohr = 0.05*(1:4096)';

Ngrid = size(g,1);

Scheme = eval(SchemeGenerator);
dR = Scheme(1).dR;

gridBohr = dR*(1:Ngrid)';
grid = gridBohr * unit2unit(DU,'Bohr',output_units.Distance);%DTab(DI('Bohr'),DI(output_units.Distance));

[mu_HNC, mu_G, mu_HNCB, mu_Nt, mu_PWC, VUA, mu_KH ] =standalone_mu_calc(gamma,c,fname,input_units,output_units,K_HK_File,['PARAMETERS;' xparam],Ngrid,dR,32768,0.00625);

Mu=struct('HNC',mu_HNC,'G',mu_G,'HNCB',mu_HNCB,'PW',mu_Nt,'PWC',mu_PWC,'VUA',VUA,'KH',mu_KH );

[f_s,pot_long] = standalone_ng_calc(fname,input_units,'PARAMETERS',Ngrid,dR);