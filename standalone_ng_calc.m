function [F_S,F_L,F_L_K]=standalone_ng_calc(fname,input_units,params,N,dR)

R = dR*(1:N)';
dk = pi/(N+1)/dR;
K = dk*(1:N)';

eval(params);

[SoluteXYZ,SoluteCharges] =load_solute2(fname,input_units);


[F_S,F_L,F_L_K]=ng_fun_calc2(SoluteCharges,Charge_o,Charge_h,Cutoff_t_ng,R,K);


