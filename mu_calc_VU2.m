function [VUA,VU,corr]=mu_calc_VU2(H_K,H_SL_K,beta,density)

[N,n]=size(H_SL_K);

H_SL0=sum(H_SL_K(1,1:end))/n;
VU=1+ density*(H_K(1,1)-H_SL0); 
VUA=VU / density*(0.52918^3);


a=-5.03;
corr=a*VU/ beta*627.50956;