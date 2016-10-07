function  [mu_HNC, mu_G, mu_HNCB, mu_Nt, mu_PWC, VUA, mu_KH ] =standalone_mu_calc(gamma,c_s,fname,input_units,output_units,KHK_file,params,N0,dR0,N,dR)



eval(params);

[SoluteXYZ,SoluteCharges,SoluteSigma,SoluteEpsilon] =load_solute2(fname,input_units);

[IN_W,DISTANCES] = build_distance_index2(SoluteXYZ,Tape,Level);


KHK = load(KHK_file);

if isstruct(KHK)
    KHK = struct2array(KHK);
end

K0 = KHK(:,1);
HK0 = KHK(:,2:end);

[R,K,H_K, W_SL_SL, ...
 F_S,F_L,F_L_K,LJ_SL_SV,br,EXPONENTA,u_short...
 chi_k_oo,chi_k_oh,chi_k_hh,chi_k_hh_i,chi_k_hh_all ...
] =...
prepare_Grid_Dependent(N0,dR0,DISTANCES,K0,HK0,...
                       SoluteCharges,SoluteSigma,SoluteEpsilon, ...
                        r_oh,r_hh,...
                       Sigma_o,Epsilon_o,Charge_o,Sigma_h,Epsilon_h,Charge_h, ...
                       Bulk_beta,density,...
                       LJ_SL_SV_Bridge_Mix,Cutoff_t_ng,Cutoff_big_pot,vdw_rules,potential_transformator); 

save chi_k.sav chi_k_oo chi_k_oh chi_k_hh_all


 [mu_HNC, mu_G, mu_HNCB, mu_Nt, mu_PWC, VUA,mu_KH ] = mu_calc3(KHK_file, gamma,c_s,SoluteCharges,Charge_o,Charge_h,Cutoff_t_ng,...
     dR0,N,dR,br,Bulk_beta,density,r_hh,r_oh);

 %[DI,EI,DTab,ETab] = build_unit_index; 
 [DU,EU] = unit_names;
 
 kcalMol_to_EU = unit2unit(EU,'kcal/mol',output_units.Energy);%ETab(EI('kcal/mol'),EI(output_units.Energy));
 
 mu_HNC = mu_HNC * kcalMol_to_EU;
 mu_G = mu_G * kcalMol_to_EU;
 mu_HNCB = mu_HNCB *kcalMol_to_EU;
 mu_Nt = mu_Nt * kcalMol_to_EU;
 mu_PWC = mu_PWC * kcalMol_to_EU;
 mu_KH = mu_KH * kcalMol_to_EU;
 
 Angstr_to_DU = unit2unit(DU,'Angstr',output_units.Distance);
 
 VUA = VUA * Angstr_to_DU^3;
