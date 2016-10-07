function mu_G = mu_calc_G2(H_SL,C_L,R,beta,density)

[N,n]=size(H_SL);

NAtom = n/2;
dR = R(2)-R(1);
r = R*ones(1,NAtom);

koeff = 4*pi*density*dR / beta*627.50956;

mu_og=(-C_L(:,1:2:end)-0.5*H_SL(:,1:2:end).*C_L(:,1:2:end));
mu_hg=(-C_L(:,2:2:end)-0.5*H_SL(:,2:2:end).*C_L(:,2:2:end));
mu_tg=mu_og+2*mu_hg;

mu_allgi = mu_tg.*r.^2;
mu_Gi = koeff * sum( mu_allgi);
mu_G = sum(mu_Gi);