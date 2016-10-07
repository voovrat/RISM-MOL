function mu_KH = mu_calc_KH(H_SL,C_L,R,beta, density )

[N,n]=size(H_SL);

NAtom = n/2;

ThH_SL = H_SL;
ThH_SL(H_SL>0) = 0;

dR = R(2)-R(1);
r = R*ones(1,NAtom);

koeff = 4*pi*density*dR / beta*627.50956;

mu_o=(0.5*ThH_SL(:,1:2:end).^2-C_L(:,1:2:end)-0.5*H_SL(:,1:2:end).*C_L(:,1:2:end));
mu_h=(0.5*ThH_SL(:,2:2:end).^2-C_L(:,2:2:end)-0.5*H_SL(:,2:2:end).*C_L(:,2:2:end));

mu_t=mu_o+2*mu_h;

mu_alli = mu_t.*r.^2;
mu_KHi = koeff * sum(mu_alli);
mu_KH = sum(mu_KHi);
