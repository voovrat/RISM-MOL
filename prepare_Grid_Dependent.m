function [R,K,H_K, W_SL_SL, ...
         F_S,F_L,F_L_K,LJ_SL_SV,br,EXPONENTA,u_short,...
        chi_k_oo,chi_k_oh,chi_k_hh,chi_k_hh_i,chi_k_hh_all ...
          ] =...
prepare_Grid_Dependent(N,dR,DISTANCES,K0,HK0,...
                       SoluteCharges,SoluteSigma,SoluteEpsilon, ...
                       r_oh,r_hh, ...
                       Sigma_o,Epsilon_o,Charge_o,Sigma_h,Epsilon_h,Charge_h, ...
                       Bulk_beta,density,...
                       LJ_SL_SV_Bridge_Mix,Cutoff_t_ng,Cutoff_big_pot,vdw_rules,potential_transformator)
                   
R = dR*(1:N)';
dk = pi/(N+1)/dR;
K = dk*(1:N)';
H_K = change_grid2(HK0,K0,K,HK0(1,:),zeros(1,size(HK0,2)));

W_SL_SL = build_W_SL_SL(K,DISTANCES);

W_oh = sin(K*r_oh)./(K*r_oh);
W_hh = sin(K*r_hh)./(K*r_hh);

[chi_k_oo,chi_k_oh,chi_k_hh,chi_k_hh_i,chi_k_hh_all] = build_Chi(H_K,W_oh,W_hh,density);

[F_S,F_L,F_L_K]=ng_fun_calc2(SoluteCharges,Charge_o,Charge_h,Cutoff_t_ng,R,K);

[LJ_SL_SV,br]=build_LJ_SL_SV(SoluteSigma,SoluteEpsilon,Sigma_o,Epsilon_o,Sigma_h,Epsilon_h,R,W_oh,W_hh,Bulk_beta,LJ_SL_SV_Bridge_Mix,vdw_rules);
             
[EXPONENTA,u_short] = short_range_pot2(LJ_SL_SV,F_S,Bulk_beta,Cutoff_big_pot);


 if strcmp(potential_transformator,'XMSA')
        
       %U_S = u_short ;
       
       M=size(u_short,2);
       for j=1:M
           i_max=min(find(EXPONENTA(:,j)>1));
          %[~,i_max]=max(EXPONENTA(:,j));
          %i_max=round(i_max*2);
          %r_max = R(i_max);
           %u_short(i_max+1:end) = 0;
           
           U_L = u_short(:,j);
           U_L(1:i_max)=0;
           u_short(i_max+1:end,j)=0;
                      
           
           F_L_K(:,j) = F_L_K(:,j) + d3fft2(U_L/Bulk_beta,R);
           F_L(:,j) = F_L(:,j)+U_L/Bulk_beta;
           F_S(:,j) = F_S(:,j)-U_L/Bulk_beta;
       end
       
       EXPONENTA = exp(-u_short);
 else   
        eval(potential_transformator);   
 end
                   

                   