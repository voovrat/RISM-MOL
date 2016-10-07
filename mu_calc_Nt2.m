function mu_Nt = mu_calc_Nt2(H_K, H_SL_K, C_L,beta,density,R,r_hh,r_oh)

[N,n]=size(H_SL_K);

dR = R(2) - R(1);

    dk = pi/(N+1)/dR;
    k = dk*(1:N)';

NAtom = n/2;

r = R*ones(1,NAtom);

koeff = 4*pi*density*dR / beta*627.50956;

sinrk_hh = sin(r_hh*k)./(r_hh*k);
sinrk_oh = sin(r_oh*k)./(r_oh*k);

DD = 1 + sinrk_hh - 2*sinrk_oh.^2;

HKDD_OO = H_K(:,1)./DD;
HKDD_OH = H_K(:,2)./DD;
HKDD_HH = H_K(:,3)./DD;

DelOO = HKDD_OO.*(1 + sinrk_hh) - 2*HKDD_OH.*sinrk_oh ;
DelOH = HKDD_OH.*(1 + sinrk_hh) - 2*HKDD_HH.*sinrk_oh ;
DelHO = HKDD_OH - HKDD_OO.*sinrk_oh;
DelHH = -HKDD_OH.*sinrk_oh + HKDD_HH;

DelOOn = DelOO*ones(1,NAtom);
DelOHn = DelOH*ones(1,NAtom);
DelHOn = DelHO*ones(1,NAtom);
DelHHn = DelHH*ones(1,NAtom);

D = density;

Del_KOO =   D * H_SL_K(:,1:2:end).*DelOOn;
Del_KOH = 2*D * H_SL_K(:,1:2:end).*DelOHn;
Del_KHO  =  D * H_SL_K(:,2:2:end).*DelHOn;
Del_KHH  =2*D * H_SL_K(:,2:2:end).*DelHHn;

Del_ROO = d3ifft2(Del_KOO,R);
Del_ROH = d3ifft2(Del_KOH,R);
Del_RHO = d3ifft2(Del_KHO,R);
Del_RHH = d3ifft2(Del_KHH,R);

CR=zeros(N,n);

CR(:,1:2:end) =    C_L(:,1:2:end).*Del_ROO + C_L(:,2:2:end).*Del_ROH;
CR(:,2:2:end) = 2*( C_L(:,2:2:end).*Del_RHH + C_L(:,1:2:end).*Del_RHO) ;

CRT = CR(:,1:2:end) + CR(:,2:2:end);

C_T = C_L(:,1:2:end) + 2 * C_L(:,2:2:end);

mu_alln = (-C_T + 0.5*CRT).*r.^2;
mu_N = koeff * sum(mu_alln);
mu_Nt = sum(mu_N);
