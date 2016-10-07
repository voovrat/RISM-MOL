function LJ_SL_SV = lj_potential_oplsaa2(C_o6, C_o12, C_h6, C_h12, R)

L = length(C_o6);
N =length(R);

R6 = R.^6;
R12 = R.^12;

r6  = R6*ones(1,L);
r12  =R12*ones(1,L);

C_o6N = ones(N,1)*C_o6';
C_o12N = ones(N,1)*C_o12';


C_h6N = ones(N,1)*C_h6';
C_h12N = ones(N,1)*C_h12';

LJ_SL_SV = zeros(N,2*L);

LJ_SL_SV(:,1:2:end) = C_o6N./r6 + C_o12N./r12;
LJ_SL_SV(:,2:2:end) = C_h6N./r6 + C_h12N./r12;


for j=1:L
    
    I = find( isinf(LJ_SL_SV(:,j)) | isnan(LJ_SL_SV(:,j)));
   
    if ~isempty(I)
        M = max(I);
        LJ_SL_SV(1:M,j) = LJ_SL_SV(M+1,j);
    end
    
end

