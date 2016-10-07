function [mu_HNC, mu_G, mu_HNCB, mu_Nt, mu_PWC, VUA, mu_KH ] = mu_calc3(KHK_file, Gamma0,c_s0,charge_sl,charge_o,charge_h,t_cut,dR0,N1,dR1,br,beta,density,r_hh,r_oh,potential_transformator)


[N0,n] = size(Gamma0);

    dk1 = pi/(N1+1)/dR1;
    k1 = dk1*(1:N1)';


R0 = dR0*(1:N0)';
R1 = dR1*(1:N1)';


K_HK = load(KHK_file);

if(isstruct(K_HK))
    K_HK = struct2array(K_HK);
end

k0 = K_HK(:,1);
H_K0 = K_HK(:,2:end);

H_K = change_grid2(H_K0,k0,k1,H_K0(1,:),zeros(1,3),k0(2)-k0(1));

Gamma1 = change_grid2(Gamma0, R0, R1, Gamma0(1,:) , zeros(1,n) , dR0);
c_s1 = change_grid2(c_s0, R0,R1, c_s0(1,:), zeros(1,n), dR0 );

br1 = change_grid2(br,R0,R1,br(1,:),zeros(1,n),dR0);

[F_S, F_L ] = ng_fun_calc2(charge_sl, charge_o, charge_h, t_cut, R1, k1);


C_L  = c_s1 - beta*F_L;
H_SL=Gamma1+c_s1; 

mu_HNC = mu_calc_HNC2(H_SL,C_L,R1,beta,density);

mu_G = mu_calc_G2(H_SL,C_L,R1,beta,density);

mu_HNCB = mu_calc_HNCB2(H_SL,C_L,R1,beta,density,br1);

H_SL_K = d3fft2( H_SL,R1);

mu_Nt = mu_calc_Nt2(H_K, H_SL_K, C_L,beta,density,R1,r_hh,r_oh);

[VUA,VU,corr]=mu_calc_VU2(H_K,H_SL_K,beta,density);

mu_PWC=mu_Nt+corr;

mu_KH = mu_calc_KH(H_SL,C_L,R1,beta,density);
