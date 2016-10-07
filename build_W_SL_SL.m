function W_SL_SL = build_W_SL_SL(k,DISTANCES)

L = length(DISTANCES);
N = length(k);

K = k*ones(1,L-1);

if L==1
    W_SL_SL = ones(N,1);
    return
end

DN = ones(N,1)*DISTANCES(2:end)';

W_SL_SL = [ ones(N,1) sin(DN .* K )./( DN.*K) ]; 
