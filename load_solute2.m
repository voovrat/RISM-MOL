function [SoluteXYZ,SoluteCharges,SoluteSigma,SoluteEpsilon] =load_solute2(filename,units)
% so far it works onle for OPLSAA force-field
% energy_units -  units for epslion: 'kJ' - kJ/mol; 'kcal' - kcal/mol
%eval(['load  -ascii ', filename, '.txt']);

%[DI,EI,DTab,ETab] = build_unit_index;
[DU,EU]=unit_names;

System = load(filename);
if isstruct(System)
    System = struct2array(System);
end

SoluteXYZ = System(:,1:3) * unit2unit(DU,units.Distance,'Bohr');%DTab(DI(units.Distance),DI('Bohr'));
SoluteSigma = System(:,4) *  unit2unit(DU,units.Distance,'Bohr');%DTab(DI(units.Distance),DI('Bohr'));
SoluteEpsilon = System(:,5) * unit2unit(EU,units.Energy,'Hartree'); %ETab(EI(units.Energy),EI('Hartree'));
SoluteCharges=System(:,6);


sigma_h_single=(0.4/0.5291772);  % LJ parameters for a single H with standard zero LJ parameters; 0.4 Angstroem
 epsilon_h_single=(0.2*300*3.1668e-006); %  LJ parameters a single H with standard zero LJ parameters; 0.2 kBT(300K);
%!!! Only for the test
 
I = find((SoluteSigma == 0) & (SoluteEpsilon==0));

if ~isempty(I)
    warning(' Zero Hydrogen Parameters!!!');
    SoluteSigma(I) = sigma_h_single;
    SoluteEpsilon(I) = epsilon_h_single;
end

% kJ_mole_to_hartree=3.808884091663899e-04; % 1 kJ/mole = 3.808884091663899*10^(-04) hartree 
% kC_mole_to_kJ= 0.239005736; % 1 kCal/mol = 0.239005736 kJ/mole;
% bohr_radius=0.5291772; % Bohr radius in Angstroems;
% %eval(['System=',filename,';']);
% [SolSize,j]=size(System);
% 
% if j~=6
%     warning('Potential parameters are missed. Probably, only coordinates are contained in the solute data.')
% end;
% 
% SoluteXYZ = System(:,1:3)/bohr_radius;

% switch energy_units
% case 'kcal'
% SoluteEpsilon = System(:,5)/kC_mole_to_kJ*kJ_mole_to_hartree;
%     otherwise % it means kJ/mol
%      SoluteEpsilon = System(:,5)*kJ_mole_to_hartree;
% end

% 
% SoluteCharges=System(:,6);




