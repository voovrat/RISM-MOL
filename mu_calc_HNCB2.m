function mu_HNCB = mu_calc_HNCB2(H_SL,C_L,R,beta,density,br)

[N,n]=size(H_SL);

NAtom = n/2;

dR = R(2)-R(1);
r = R*ones(1,NAtom);

koeff = 4*pi*density*dR / beta*627.50956;

mu_ob=(0.5*H_SL(:,1:2:end).^2-C_L(:,1:2:end)-0.5*H_SL(:,1:2:end).*C_L(:,1:2:end)+(H_SL(:,1:2:end)+1).*(exp(-abs(br(:,1:2:end)))-1));
mu_hb=(0.5*H_SL(:,2:2:end).^2-C_L(:,2:2:end)-0.5*H_SL(:,2:2:end).*C_L(:,2:2:end)+(H_SL(:,2:2:end)+1).*(exp(-abs(br(:,2:2:end)))-1));
mu_tb=mu_ob+2*mu_hb;



% for j_mu_calc=1:NAtom; for i_mu_calc=1:Bulk.Grid_N;  mu_allbi(i_mu_calc,j_mu_calc)=mu_tb(i_mu_calc,j_mu_calc)*R(i_mu_calc)^2;end;
%     mu_HNBi(j_mu_calc)=4*pi*Bulk.density*sum(mu_allbi(1:end,j_mu_calc))*Bulk.Grid_dR/Bulk.beta*627.50956;
% end;
% 
% mu_HNCB=sum(mu_HNBi);

mu_allbi = mu_tb.*r.^2;
mu_HNBi = koeff* sum(mu_allbi);
mu_HNCB = sum(mu_HNBi);