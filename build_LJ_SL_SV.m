function [LJ_SL_SV,br] = build_LJ_SL_SV(SoluteSigma,SoluteEpsilon,Sigma_o,Epsilon_o,Sigma_h,Epsilon_h,R,W_OH,W_HH,Bulk_beta,K_br,vdw_rules) 

[sigma_new_o, epsilon_new_o, sigma_new_h, epsilon_new_h] = feval(vdw_rules,SoluteSigma, SoluteEpsilon, Sigma_o, Epsilon_o, Sigma_h, Epsilon_h);

[C_o6,C_o12,C_h6,C_h12] = SigmaEpsilon2C6C12(sigma_new_o, epsilon_new_o, sigma_new_h, epsilon_new_h);

LJ_SL_SV0 = lj_potential_oplsaa2(C_o6,C_o12,C_h6,C_h12, R);

 
    br = buildRepulsiveBridge(C_o12,C_h12,R,W_OH,W_HH,Bulk_beta);
  
    LJ_SL_SV = LJ_SL_SV0 + K_br*br/Bulk_beta;



