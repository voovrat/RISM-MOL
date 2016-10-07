function [gamma_out,c_out,Err,timing,LJ_SL_SV,EXPONENTA,u_short,F_S,F_L,F_L_K]=run_NMGM_iterations(filename,input_units,KHK_file,params,Scheme,R0,gamma0)



callback = 0;
DownCallback=0;
UpCallback=0;
eval(params);

[SoluteXYZ,SoluteCharges,SoluteSigma,SoluteEpsilon] =load_solute2(filename,input_units);

[IN_W,DISTANCES] = build_distance_index2(SoluteXYZ,Tape,Level);

KHK = load(KHK_file);

if isstruct(KHK)
    KHK = struct2array(KHK);
end

K0 = KHK(:,1);
HK0 = KHK(:,2:end);

S =  struct('N',NaN,'dR',NaN,'Niter_pre',NaN,'Niter_post',NaN,'eps',NaN,'mu',NaN);
SchemeCell = struct2cell(orderfields(Scheme,S))';

Nc = SchemeCell(:,1);
dRc = SchemeCell(:,2);

Nn = cell2mat(Nc);
dRn = cell2mat(dRc);


Ng_function

 [Rn,Kn,H_Kn, W_SL_SLn, ...
         F_Sn,F_Ln,F_L_Kn,LJ_SL_SVn,brn,EXPONENTAn,u_shortn...
        chi_k_oon,chi_k_ohn,chi_k_hhn,chi_k_hh_in,chi_k_hh_alln ...
          ] =...
prepare_MGM_Grid_Dependent(Nn,dRn,DISTANCES,K0,HK0,...
                       SoluteCharges,SoluteSigma,SoluteEpsilon, ...
                       r_oh,r_hh, ...
                       Sigma_o,Epsilon_o,Charge_o,Sigma_h,Epsilon_h,Charge_h, ...
                       Bulk_beta,density,...
                       LJ_SL_SV_Bridge_Mix,Cutoff_t_ng,Cutoff_big_pot,vdw_rules,...
                       potential_transformator,Ng_function,Short_Range_FN,MSACutFactor);

             
L = size(SoluteXYZ,1);                   

if ~ismember('NMGM_StartLevel',who) || NMGM_StartLevel == 0
  NMGM_StartLevel = length(Scheme);
end

MaxLev = NMGM_StartLevel;
  
if nargin<7 || sum(sum(isnan(R0) ))
    gamma_out = zeros(Nn(MaxLev),2*L);                   
else
    gamma_out = change_grid2(gamma0,R0,Rn{MaxLev},gamma0(1,:),zeros(1,2*L));
end



t0 = mytic;
timing = struct;

if ismember('Write_Protocol_dGamma',who) && Write_Protocol_dGamma
    timing.write_dGamma = 1;  
else
       timing.write_dGamma = 0; 
end

if ~ismember('Lambda_CoarseSteps',who) 
    Lambda_CoarseSteps=0;
end

   timing.Lambda = Lambda_CoarseSteps;
   timing.A = CoarseSteps_Treshold;
    timing.Nested = Nested_Iterations;

C = struct2cell(Scheme)';
FirstLev = find(cell2mat(C(:,end))>1000);

if isempty(FirstLev)
    FirstLev=1;
else
    FirstLev = FirstLev(end);
end

%Scheme(FirstLev).Niter_post = 1;


for lev = MaxLev:-1:FirstLev
   
   
  %  for MGMit = 1:nit
    [gamma_out,c_out,Err,timing] = do_NMGM_iteration(Scheme,lev,gamma_out,0,Rn,Kn, IN_W, EXPONENTAn, u_shortn, F_L_Kn, ...
        Bulk_beta,W_SL_SLn, chi_k_oon, chi_k_ohn,chi_k_hh_alln, Lambda_Mix, closure, ...
        callback,DownCallback,UpCallback,RCorr,timing);
  
    
    if sum(sum(isnan(gamma_out)) )
        return
    end
%    end

  
   if lev>1
        RC = cell2mat(Rn(lev));
        RF = cell2mat(Rn(lev-1));

        gamma_old = gamma_out;
        gamma_out = change_grid2(gamma_out,RC,RF,gamma_out(1,:),zeros(1,2*L));
   end
   
end

%%% prolongation

for lev = FirstLev-1:-1:1
    
  
    Scheme1=Scheme;
    
    Dz = FirstLev-lev;
    
    Scheme1 = [Scheme(Dz+1:end); Scheme(1:Dz)  ];
    
    for i=1:length(Scheme1)
        Scheme1(i).N=Scheme1(i).N * 2^(FirstLev-lev);
    end
   
   SchemeCell = struct2cell(orderfields(Scheme1,S))';

Nc1 = SchemeCell(:,1);
dRc1 = SchemeCell(:,2);

Nn1 = cell2mat(Nc1);
dRn1 = cell2mat(dRc1);

 [Rn1,Kn1,H_Kn1, W_SL_SLn1, ...
         F_Sn1,F_Ln1,F_L_Kn1,LJ_SL_SVn1,brn,EXPONENTAn1,u_shortn1...
        chi_k_oon1,chi_k_ohn1,chi_k_hhn1,chi_k_hh_in1,chi_k_hh_alln1 ...
          ] =...
prepare_MGM_Grid_Dependent(Nn1,dRn1,DISTANCES,K0,HK0,...
                       SoluteCharges,SoluteSigma,SoluteEpsilon, ...
                       r_oh,r_hh, ...
                       Sigma_o,Epsilon_o,Charge_o,Sigma_h,Epsilon_h,Charge_h, ...
                       Bulk_beta,density,...
                       LJ_SL_SV_Bridge_Mix,Cutoff_t_ng,Cutoff_big_pot,vdw_rules,potential_transformator,Ng_function,Short_Range_FN);
           

   if   isfield(timing,'prev_dGamma')              
            timing=rmfield(timing,'prev_dGamma');
   end
   if isfield(timing,'nit');
       timing = rmfield(timing,'nit');
   end
   
   
   if isfield(timing,'prev_dGammaErr')
            timing=rmfield(timing,'prev_dGammaErr');
            Scheme1(end-Dz).Niter_pre = timing.CoarsestIter;
   end
    
   [gamma_out,c_out,Err,timing] = do_NMGM_iteration(Scheme1,lev,gamma_out,0,Rn1,Kn1, IN_W, EXPONENTAn1, u_shortn1, F_L_Kn1, ...
        Bulk_beta,W_SL_SLn1, chi_k_oon1, chi_k_ohn1,chi_k_hh_alln1, Lambda_Mix, closure, ...
        callback,DownCallback,UpCallback,RCorr,timing);
 
    
 %   W_SL_SL = cell2mat(W_SL_SLn(lev)) ;

 %   EXPONENTA = cell2mat(EXPONENTAn(lev));
 %   u_short = cell2mat(u_shortn(lev));
 %   F_L_K = cell2mat(F_L_Kn(lev));
    
 %   chi_k_oo = cell2mat(chi_k_oon(lev));
 %   chi_k_oh = cell2mat(chi_k_ohn(lev));
 %   chi_k_hh_all=cell2mat(chi_k_hh_alln(lev));
    
 %   R=Rn{lev};
 %   K=Kn{lev};
   
 %   Niter = Scheme(lev).Niter_post;
 %   eps = Scheme(lev).eps;
     
 %   t0 = tic ;
      
    
    
 %   [gamma_out,c_out,niter,Err] = CALC_LONGS3(R,K,gamma_out,0,IN_W,EXPONENTA, ...
 %                u_short,F_L_K,Bulk_beta,W_SL_SL,chi_k_oo,chi_k_oh,chi_k_hh_all,Lambda_Mix,...
 %            closure,Niter,eps,callback,RCorr);

 %   dt = toc(t0);
  %  timing.t_calc(lev) = timing.t_calc(lev) + dt;
  %  timing.protocol= [ timing.protocol; struct('lev',lev,'niter',niter,'Err',Err,'time',dt) ];   
         
   if lev>1
        RC = cell2mat(Rn(lev));
        RF = cell2mat(Rn(lev-1));

        gamma_old = gamma_out;
        gamma_out = change_grid2(gamma_out,RC,RF,gamma_out(1,:),zeros(1,2*L));
   end
     
end

timing.t_total_calc = mytoc(t0);


LJ_SL_SV = LJ_SL_SVn{1};
EXPONENTA = EXPONENTAn{1};
u_short = u_shortn{1};
F_S = F_Sn{1};
F_L = F_Ln{1};
F_L_K = F_L_Kn{1};