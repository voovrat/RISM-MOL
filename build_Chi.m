function [chi_k_oo,chi_k_oh,chi_k_hh,chi_k_hh_i,chi_k_hh_all] = build_Chi(HK,W_oh,W_hh,density)
    
chi_k_oo=1+ HK(:,1)*density;
chi_k_oh=W_oh + HK(:,2)*density;
chi_k_hh=W_hh + HK(:,3)*density;% internal function H1H2 or H2H1;
chi_k_hh_i= 1 + HK(:,3)*density; % internal function H1H1 or H2H2;
chi_k_hh_all=chi_k_hh + chi_k_hh_i; % internal function H1H2 AND H1H1;