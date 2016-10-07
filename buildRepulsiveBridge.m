function br = buildRepulsiveBridge(C_o12,C_h12,R,W_OH,W_HH,Bulk_beta)
%  
%  br = buildRepulsiveBridge(C_o12,C_h12,R,W_OH,W_HH,Bulk_beta)
%   
%   repulsive hydrogen and oxygen bridge by Kovalenko
%
%    exp(-B[s,alpha]) = Prod[nu<>alpha] <W[alpha,nu] * exp(-beta*4*epsilon[s,nu](sigma[s,nu]/r)^12) >
%   
% Notes to implementation:
%   1. 4 *epsilon * (sigma/r)^12 = C12/r^12
%
%   2. exp(-B[s,o]) = < W_oh * exp(-beta*C12[s,h]/r^12) >^2
%      exp(-B[s,h]) = <W_hh * exp(-beta*C12[s,h]/r^12)>*<W_oh * exp(-beta*C12[s,o]/r^12) >   
%

L=length(C_o12);
N = length(R);

R12 = R.^12;
r12 = R12 * ones(1,L);

C_o12N = ones(N,1)*C_o12';
C_h12N = ones(N,1)*C_h12';
W_OHN = W_OH*ones(1,L);
W_HHN = W_HH*ones(1,L);


% <exp(-beta/r)*w> = <(exp(-beta/r)-1)*w>+w 

exp_h = exp(-Bulk_beta*C_h12N./r12)-1;
exp_o = exp(-Bulk_beta*C_o12N./r12)-1;

exp_h_k = d3fft2(exp_h,R);
exp_o_k = d3fft2(exp_o,R);

exp_br_o = (d3ifft2(W_OHN.*exp_h_k,R)+1).^2;
exp_br_h = (d3ifft2(W_HHN.*exp_h_k,R)+1) .* (d3ifft2(W_OHN.*exp_o_k,R)+1);

%exp_br_o = (d3ifft2(W_OHN.*exp_h_k,R)+1);
%exp_br_h = (d3ifft2(W_OHN.*exp_o_k,R)+1);


br = zeros(N,2*L);

% negative values should appear only because of truncation errors
exp_br_o(exp_br_o<1e-12) = 1e-12;
exp_br_h(exp_br_h<1e-12) = 1e-12;

% at the last part of exp_br appears regions of decay
% that's because of finite domain
%exp_br_o([ false ;diff(exp_br_o)<-1e-5]) = 1;
%exp_br_h([ false ;diff(exp_br_h)<-1e-5]) = 1;
exp_br_o(N/2:end,:)=1;
exp_br_h(N/2:end,:)=1;

br_o = -log(exp_br_o);
br_h = -log(exp_br_h);

br(:,1:2:end) = br_o;
br(:,2:2:end) = br_h;