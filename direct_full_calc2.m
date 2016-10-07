function [rGamma_out, r_cs] = direct_full_calc2(r,rGamma_in, IN_W, EXPONENTA, u_short, k_fk, W_SL_SL, chi_k_oo, chi_k_oh,chi_k_hh_all, lambda, closure,NCorr)
% size(r) = N x N_SL*2
% 


N_SL = size(rGamma_in,2)/2;
N = size(r,1);
dr = r(2,1) - r(1,1);

dk = pi/(N+1)/dr;

r_cs = feval(closure,r,rGamma_in, EXPONENTA, u_short); 

I = isnan(r_cs);
r_cs(I) = -r(I) - rGamma_in(I);

k_ck = 4*pi*dst2(r_cs)*dr - k_fk; % k_ck = k_ck0 - k_fk. to avoid substraction inside the cycle...

kGamma_k = zeros(size(rGamma_in));

for sol=1:N_SL;
    for j=1:N_SL;
        
        p = IN_W(sol,j);
        
        if p
            
            kGamma_k(:,2*sol-1)=kGamma_k(:,2*sol-1)+...
                W_SL_SL(:,p).* ( k_ck(:,2*j-1).*chi_k_oo  +  2*k_ck(:,2*j) .* chi_k_oh); 
            
            kGamma_k(:,2*sol)=kGamma_k(:,2*sol)+...
                W_SL_SL(:,p).*( k_ck(:,2*j-1).*chi_k_oh + k_ck(:,2*j).*chi_k_hh_all );

            
        end
    end
end;

kGamma_k = kGamma_k - k_ck -  k_fk; % k_ck = ( cf prev comment)= k_ck0 - k_fk.
                                    % gamma = <w*(k_ck0-k_fk)*chi> -
                                    % (kc_k0-k_fk)-k_fk =
                                    % <w*(k_ck0-k_fk)*chi>  - k_ck0

r_gamma_tmp =  (1/2/pi^2*dst2(kGamma_k)*dk);

if NCorr
    r_gamma_tmp(1:NCorr,:) = r_gamma_tmp(1:NCorr,:) + r_cs(1:NCorr,:) + r(1:NCorr,:);
end
                                    
%rGamma_out = (1-lambda)*rGamma_in  + lambda*(1/2/pi^2*dst2(kGamma_k)*dk);

rGamma_out = (1-lambda)*rGamma_in  + lambda*r_gamma_tmp;