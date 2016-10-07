function [gamma_out,c_out,iter,Err] = CALC_LONGS3(R,K,gamma_in,d, IN_W, EXPONENTA, u_short, f_k,Bulk_beta,W_SL_SL, chi_k_oo, chi_k_oh,chi_k_hh_all, lambda, closure,max_iter,eps,callback,RCorr)
%
% we expect L(gamma) := G(C(gamma)) - gamma
% ( where G(C(gamma)) == direct_full_calc2(gamma) )
% and now solve L(gamma) = d
%
% so, iterations are
%
% gamma = G(C(gamma))-d
%
%
%

if nargin<18
    callback=0;
end

if nargin<19
    RCorr=0
end

c_out=0;
iter=1;

L = size(gamma_in,2);

dR = R(2)-R(1);

if RCorr
    NCorr = round(RCorr/dR);
else
    NCorr = 0;
end

r = R*ones(1,L); 

dk=K(2)-K(1);

k = K*ones(1,L);

k_fk = Bulk_beta*k.*f_k;

 rGamma_in = gamma_in.*r;
 
% rGamma_out = direct_full_calc2(r,dk,rGamma_in, IN_W, EXPONENTA, u_short, k_fk, W_SL_SL, chi_k_oo, chi_k_oh,chi_k_hh_all, lambda, closure) - r.*d;
rGamma_out =direct_full_calc2(r,rGamma_in, IN_W, EXPONENTA, u_short, k_fk, W_SL_SL, chi_k_oo, chi_k_oh,chi_k_hh_all, lambda, closure,NCorr) - r.*d;

Err = ones(10000,1)*inf;

 for iter=1:max_iter
    
  rGamma_old = rGamma_out;
     
 %[rGamma_out, rc_s] =direct_full_calc2(r,dk,rGamma_out, IN_W, EXPONENTA, u_short, k_fk, W_SL_SL, chi_k_oo, chi_k_oh,chi_k_hh_all, lambda, closure);
 [rGamma_out, rc_s] =direct_full_calc2(r,rGamma_out, IN_W, EXPONENTA, u_short, k_fk, W_SL_SL, chi_k_oo, chi_k_oh,chi_k_hh_all, lambda, closure,NCorr);
 
 
 if sum(sum(isnan(rGamma_out)))
     Err= NaN;
     break
 end
    
 
 rGamma_out = rGamma_out - r.*d;
 
 %Err = sqrt(sum(sum((rGamma_out - rGamma_old).^2)))*dR/L;

 if isnan(iter)
     db=1;
 else
     
 Err(iter) =  count_err(rGamma_out./r ,rGamma_old./r,dR);%mean( sqrt(sum(rGamma_out./r - rGamma_old./r).^2)*dR );
 end
 
 if callback
     eval(callback);
 end
    
 if Err(iter)<eps
        break;
 end

 end

 gamma_out = rGamma_out./r;
if max_iter
 c_out = rc_s./r;
end

if isempty(iter)
    iter=0; 
end

if ~isnan(Err)
  Err = Err(1:iter);
else
    db=1;
end
