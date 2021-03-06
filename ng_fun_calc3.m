function [F_S, F_L, F_L_K] = ng_fun_calc3(charge_sl, charge_o, charge_h, t, R, K)
% Calculation of functions for the Ng-renormalisation;
% CHARGE_SL - array with partial charges of the solute;
% CHARGE_O - Oxygen's partial charge (with sighn!);
% CHARGE_H - Hydrogen's partial charge;
% T - mixting parameter (0.8 - 1.5 usually);
% F_S 'Short-range Function'
% F_L 'Long-range Function'
% F_L_K -  'Analytic Fourier Transform of F_L';
% R - Points of grid in real space;
% K - Points of grid in reciprocal space;
L = length(charge_sl);
N = length(R);

z = zeros(1,2*L);

z(1:2:end) = charge_o*charge_sl';
z(2:2:end) = charge_h*charge_sl';

if length(t)==1

erfR=erf(t*R);
erfR1=(1-erfR);

erfRN = erfR*ones(1,2*L);
erfR1N = erfR1*ones(1,2*L);

else

    % just for debug
  erfR0=erf(t(1)*R);
erfR10=(1-erfR0);

erfRN0 = erfR0*ones(1,2*L);
erfR1N0 = erfR10*ones(1,2*L);  

% ---- end debug
    
erfR_o=erf(t(1)*R);
erfR1_o=(1-erfR_o);   

erfR_h=erf(t(2)*R);
erfR1_h=(1-erfR_h);   

erfRN_o = erfR_o*ones(1,L);
erfR1N_o = erfR1_o*ones(1,L);


erfRN_h = erfR_h*ones(1,L);
erfR1N_h = erfR1_h*ones(1,L);

erfRN = zeros(N,2*L);
erfR1N = zeros(N,2*L);

erfRN(:,1:2:end) = erfRN_o;
erfRN(:,2:2:end) = erfRN_h;

erfR1N(:,1:2:end) = erfR1N_o;
erfR1N(:,2:2:end) = erfR1N_h;


end


zN = ones(N,1)*z;

r = R*ones(1,2*L);
k2 = K.^2*ones(1,2*L);


F_S = zN.*erfR1N./r;   % 'Short Function'
%F_S = zeros(size(erfR1N));
F_L = zN.*erfRN./r;  % 'Long Function'


%N_plus = 1048576;
%R_plus = (1:N_plus)'*0.05;
%r_plus = R_plus*ones(1,2*L);
%zN_plus=ones(N_plus,1)*z;

%erfR_plus = erf(t(1)*R_plus);
%erfR1_plus=(1-erfR_plus);

%erfRN_plus = erfR_plus*ones(1,2*L);
%erfR1N_plus = erfR1_plus*ones(1,2*L);
 
%F_L_plus = zN_plus.*erfRN_plus./r_plus;

%F_L_K_d_plus=d3fft2(F_L_plus,R_plus);

%plot(F_L_K_d_plus(:,1)),xlim([0 100])






if length(t)==1
    
F_L_K = 4*pi*zN.*exp(-k2/(4*t^2))./k2; % 'Analytic Fourier Transform of F_L';

else
  
    % just for debug ---
F_L0 = zN.*erfRN0./r;
F_S0 = zN.*erfR1N0./r;
 % --- d debug
    
 F_L_K0  =  4*pi*zN.*exp(-k2/(4*t(1)^2))./k2;
    
 F_L_K_o = 4*pi*zN(:,1:2:end).*exp(-k2(:,1:2:end)/(4*t(1)^2))./k2(:,1:2:end); % 'Analytic Fourier Transform of F_L';   
 F_L_K_h = 4*pi*zN(:,2:2:end).*exp(-k2(:,2:2:end)/(4*t(2)^2))./k2(:,2:2:end); % 'Analytic Fourier Transform of F_L';   
 
 F_L_K = zeros(N,2*L);
 F_L_K(:,1:2:end) = F_L_K_o;
 F_L_K(:,2:2:end) = F_L_K_h;
 
end


Sigma12_h = (0.4+3.16)/0.52/2;
Rep_h=7.3306e-005*(Sigma12_h./R).^12;
Rep_h(Rep_h>1000) = 1000;
Rep_h_K = d3fft2(Rep_h,R);

F_tot_K_h = F_L_K(:,2) + Rep_h_K;

F_tot_h = d3ifft2(F_tot_K_h,R);
save F_tot_h.dat F_tot_h -ASCII -double
%figure,plot(R,F_tot_h),xlim([0 10]),ylim([-1 1]);