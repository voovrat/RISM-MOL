function [gamma_out,c_out,Error,timing] = do_NMGM_iteration(Scheme,lev, gamma_in,f, ...
                    Rn,Kn, IN_W, EXPONENTAn, u_shortn, F_L_Kn,Bulk_beta, ...
                W_SL_SLn, chi_k_oon, chi_k_ohn,chi_k_hh_alln, lambda, closure,...
                callback,DownCallback,UpCallback,RCorr,timing)
% L(gamma) = G(C(gamma)) - gamma = direct_full_calc(gamma) - gamma 
% we solve L(gamma) = f
% 
% NMGM(L,gamma_fine,f):
%     gamma_fine = DoIterations(L,gamma_fine,f,Niter_pre)
%     if L==max_level return gamma_fine
%
%     gamma_coarse = change_grid(gamma_fine)
%     defect_fine = G(C(gamma_fine)) - f  = direct_full_calc(gamma_fine) -gamma_fine
%     defect_coarse = change_grid(defect_fine)
%
%     d := L(gamma_coarse) - defect_coarse = G(C(gamma_coarse)) -
%     gamma_coarse - defect_coarse
%
%     v = gamma_coarse;
%
%     for i = 1:mu do NMGM(L+1, v, d )
%   
%     gamma = gamma_fine + change_grid ( v - gamma_coarse)
%
%     gamma_out = DoIterations(L,gamma,f,Niter_post); 


if nargin<22
    timing = struct;
end
   
if ~isfield(timing,'t_calc') || length(timing.t_calc)<lev
    timing.t_calc(lev) = 0;
end

if  ~isfield(timing,'protocol') 
    timing.protocol = [];
end

if ~isfield(timing,'nit')
    timing.nit=0;
end



if sum(sum(isnan(gamma_in)))
    return
end

%--- fine grid parameters
    R =cell2mat(Rn(lev));
    K = cell2mat(Kn(lev));

%    L = size(gamma_in,2);
% 
%     k = K*ones(1,L);
%     r = R*ones(1,L);
%     dk = K(2) - K(1);

    W_SL_SL = cell2mat(W_SL_SLn(lev)) ;

    EXPONENTA = cell2mat(EXPONENTAn(lev));
    u_short = cell2mat(u_shortn(lev));
    F_L_K = cell2mat(F_L_Kn(lev));

    chi_k_oo = cell2mat(chi_k_oon(lev));
    chi_k_oh = cell2mat(chi_k_ohn(lev));
    chi_k_hh_all=cell2mat(chi_k_hh_alln(lev));
%---
    
Niter_pre = Scheme(lev).Niter_pre;
Niter_post = Scheme(lev).Niter_post;
eps = Scheme(lev).eps;
mu = Scheme(lev).mu;

L=size(gamma_in,2);

t0 = mytic ;


%if length(Scheme)>1 && lev == length(Scheme)
%    eps1 = 0;
%else
%    eps1 = eps;
%end


if Scheme(lev).mu==0

    A = timing.A; 
    Lam = timing.Lambda;
    
    if Lam && timing.nit 
       eps_dGamma  = count_err(timing.prev_dGamma,f,R(2)-R(1));
      
       
       
       nstep = timing.CoarseSteps;
       delta = (timing.e1/timing.e0);
       en = timing.e0 * delta^nstep;
       
       B = en/eps_dGamma;
       
       nadd = -log(B/A)/log(delta);
       
       nopt = nstep + nadd;
        
       
       Niter_pre = round((1-Lam)*nstep + Lam*nopt);
    end
    

 
    
    if timing.Lambda
        if Niter_pre <50
            Niter_pre = 50;
        end
    
    
        if Niter_pre > 2000
            Niter_pre = 2000;
        end
    end
    
    timing.prev_dGamma = f;
    timing.CoarseSteps = Niter_pre;
    timing.nit = timing.nit+1;
end




eps1=eps;

e0=0;e1=0;

if Scheme(lev).mu==0
  [gamma_1,c_out,niter,Err] = CALC_LONGS3(R,K,gamma_in,f,IN_W,EXPONENTA,u_short,F_L_K,Bulk_beta,W_SL_SL,chi_k_oo,chi_k_oh,chi_k_hh_all,lambda,closure,1,eps1,callback,RCorr);
  [gamma_out,c_out,niter,Err] = CALC_LONGS3(R,K,gamma_1,f,IN_W,EXPONENTA,u_short,F_L_K,Bulk_beta,W_SL_SL,chi_k_oo,chi_k_oh,chi_k_hh_all,lambda,closure,Niter_pre-1,eps1,callback,RCorr);

    e0 = count_err(gamma_in,gamma_out,R(2)-R(1));
    e1 = count_err(gamma_1,gamma_out,R(2)-R(1));

    timing.e0 = e0;
    timing.e1 = e1;
    
else
  [gamma_out,c_out,niter,Err] = CALC_LONGS3(R,K,gamma_in,f,IN_W,EXPONENTA,u_short,F_L_K,Bulk_beta,W_SL_SL,chi_k_oo,chi_k_oh,chi_k_hh_all,lambda,closure,Niter_pre,eps1,callback,RCorr);
end

if Niter_pre>=1
    Error=Err(end);
else 
    Error=inf;
end

dt = mytoc(t0);
timing.t_calc(lev) = timing.t_calc(lev) + dt;
%timing.n_iter(lev) = timing.n_iter(lev)+niter;
%timing.err{lev} = [timing.err{lev} ; Err];

if timing.write_dGamma
    
  
    
    
    
    if isfield(timing,'prev_dGammaErr')
        pe=timing.prev_Err;
        pge = timing.prev_dGammaErr;
    else
        pe =0;
        pge=0;
    end
    
    
     timing.protocol= [ timing.protocol; ...
                        struct('lev',lev,'niter',niter,'Err',...
                        Err,'time',dt,'dGamma',f,...
                        'gamma_in',gamma_in,'gamma_out',gamma_out,...
                        'pe',pe,'pge',pge,'e0',e0,'e1',e1) ];
else
     timing.protocol= [ timing.protocol; struct('lev',lev,'niter',niter,'Err',Err,'time',dt) ];
end
   
if isempty(Err)
    lastErr = inf;
else
    lastErr = Err(end);
end
    

if sum(sum(isnan(gamma_out)))
    Err=Err(end);
    return;
end

if DownCallback
    eval(DownCallback);
end

if lev == length(Scheme) || Scheme(lev).mu==0


    timing.E0 = Err(1);
    timing.E1 = Err(end);
    Error=Err(end);       
        
    return
end

% --- coarse grid parameters
    RC = cell2mat(Rn(lev+1));
    KC = cell2mat(Kn(lev+1));

%     rC = RC*ones(1,L);
%     kC = KC*ones(1,L);
%     dkC = KC(2)-KC(1);

    W_SL_SLC = cell2mat(W_SL_SLn(lev+1)) ;

    EXPONENTAC = cell2mat(EXPONENTAn(lev+1));
    u_shortC = cell2mat(u_shortn(lev+1));
    F_L_KC = cell2mat(F_L_Kn(lev+1));

    chi_k_ooC = cell2mat(chi_k_oon(lev+1));
    chi_k_ohC = cell2mat(chi_k_ohn(lev+1));
    chi_k_hh_allC=cell2mat(chi_k_hh_alln(lev+1));
%---

gamma_fine = gamma_out;
%rG_C_gamma_fine = direct_full_calc2(r,dk,r.*gamma_fine, IN_W, EXPONENTA, u_short, k.*F_L_K, W_SL_SL, chi_k_oo, chi_k_oh,chi_k_hh_all, lambda, closure);
%G_C_gamma_fine = rG_C_gamma_fine./r;


t0 = mytic;
[G_C_gamma_fine,ct,niter,Err] = CALC_LONGS3(R,K,gamma_out,f,IN_W,EXPONENTA,u_short,F_L_K,Bulk_beta,W_SL_SL,chi_k_oo,chi_k_oh,chi_k_hh_all,lambda,closure,0,0,0,RCorr);

dt = mytoc(t0);


Error = count_err(G_C_gamma_fine,gamma_out,R(2)-R(1));


timing.t_calc(lev) = timing.t_calc(lev) + dt;
%timing.n_iter(lev) = timing.n_iter(lev)+niter;
%timing.err{lev} = [timing.err{lev}  ; Err];
%timing.protocol= [ timing.protocol; struct('lev',lev,'niter',niter,'Err',Err,'time',dt) ];



if sum(sum(isnan(G_C_gamma_fine)))
    return
end


L_gamma_fine = G_C_gamma_fine - gamma_fine;
defect_fine = L_gamma_fine - f;

defect_coarse = change_grid2(defect_fine,R,RC,defect_fine(1,:),zeros(1,L));


gamma_coarse = change_grid2(gamma_out,R,RC,gamma_out(1,:),zeros(1,L));

if f==0
    NC=size(RC,1);
    f_coarse = zeros(NC,L);
else
    f_coarse = change_grid2(f,R,RC,f(1,:),zeros(1,L));
end
%rG_C_gamma_coarse = direct_full_calc2(rC,dkC,rC.*gamma_coarse, IN_W, EXPONENTAC, u_shortC, kC.*F_L_KC, W_SL_SLC, chi_k_ooC, chi_k_ohC,chi_k_hh_allC, lambda, closure);
%G_C_gamma_coarse = rG_C_gamma_coarse./rC;

t0 = mytic;
[G_C_gamma_coarse,ct,niter,Err] = CALC_LONGS3(RC,KC,gamma_coarse,f_coarse,IN_W,EXPONENTAC,u_shortC,F_L_KC,Bulk_beta,W_SL_SLC,chi_k_ooC,chi_k_ohC,chi_k_hh_allC,lambda,closure,0,0,0,RCorr);

dt = mytoc(t0);
timing.t_calc(lev) = timing.t_calc(lev) + dt;
%timing.n_iter(lev) = timing.n_iter(lev)+niter;
%timing.err{lev} = [timing.err{lev} ; Err];
%timing.protocol= [ timing.protocol; struct('lev',lev,'niter',niter,'Err',Err,'time',dt) ];

if sum(sum(isnan(G_C_gamma_coarse)))
    return
end

L_gamma_coarse = G_C_gamma_coarse - gamma_coarse;

d = L_gamma_coarse - defect_coarse;

%v = G_C_gamma_coarse;%gamma_coarse;
v=gamma_coarse;

%gamma_out = G_C_gamma_fine; % one iteration was done => we are not going to lost it 

for mit = 1:mu
    if timing.Nested
        break;
    end
    
    if isempty(lastErr)
        lastErr = inf; 
        break;
    end
    
    if lastErr(end) < eps
        break
     end
     
     Scheme1=Scheme;
     
  %   if mu>1
        % dGamma = timing.prev_dGamma;
      %  db=1;
   %    if isfield(timing,'prev_dGammaErr') && timing.prev_dGammaErr>0
           
    %      Lam = timing.Lambda; 
        %  Koeff0=timing.E0/timing.prev_dGammaErr;
         %  Koeff0=timing.prev_Err/timing.prev_dGammaErr;
         % Koeff0 = timing.E0
     %    A=2.6;
         
         
         
      %     Koeff =  (1-Lam) + Lam*Koeff0;
          
       %   C=struct2cell(Scheme)';
         % mu_id = cell2mat(C(:,end));
        %  
          %CoarsestIter = round(timing.last_NiterPre*Koeff);
          
         % timing.CoarsestIter = CoarsestIter;
          
         % if CoarsestIter<50
         %     CoarsestIter=50;
         % end
         % if CoarsestIter>5000
         %     CoarsestIter=5000;
         % end
          
        %  Scheme1(mu_id==0).Niter_pre = CoarsestIter;
      % end
      
    % end
     
        [v,c,lastErr,timing] = do_NMGM_iteration(Scheme1,lev+1, v, d, Rn,Kn,  IN_W, EXPONENTAn, u_shortn, F_L_Kn, ...
        Bulk_beta,W_SL_SLn, chi_k_oon, chi_k_ohn,chi_k_hh_alln, lambda, closure,...
        callback,DownCallback,UpCallback,RCorr,timing);

    if sum(sum(isnan(v)))
        return
    end
    
 
end

d_gamma_coarse = v - gamma_coarse; %G_C_gamma_coarse;



d_gamma_fine = change_grid2(d_gamma_coarse,RC,R,d_gamma_coarse(1,:),zeros(1,L));

gamma_out = gamma_out + d_gamma_fine;

t0 = mytic;
[gamma_out,c_out,niter,Err] = CALC_LONGS3(R,K,gamma_out,f,IN_W,EXPONENTA,u_short,F_L_K,Bulk_beta,W_SL_SL,chi_k_oo,chi_k_oh,chi_k_hh_all,lambda,closure,Niter_post,eps,callback,RCorr);

if Niter_post>0
    Error=Err(end);
end

dt = mytoc(t0);
timing.t_calc(lev)=timing.t_calc(lev)+dt;
%timing.n_iter(lev) = timing.n_iter(lev)+niter;
%timing.err{lev} = [timing.err{lev} ; Err];


if timing.write_dGamma

    
         timing.protocol= [ timing.protocol; ...
                        struct('lev',lev,'niter',niter,'Err',...
                        Err,'time',dt,'dGamma',f,...
                        'gamma_in',gamma_in,'gamma_out',gamma_out,...
                        'pe',pe,'pge',pge,'e0',e0,'e1',e1) ];
    
    %    timing.protocol= [ timing.protocol; struct('lev',lev,'niter',niter,'Err',Err,'time',dt,'dGamma',f,'pe',pe,'pge',pge) ];

    
    %timing.protocol= [ timing.protocol; struct('lev',lev,'niter',niter,'Err',Err,'time',dt,'dGamma',f) ];
else
     timing.protocol= [ timing.protocol; struct('lev',lev,'niter',niter,'Err',Err,'time',dt) ];
end

if UpCallback
    eval(UpCallback);
end



