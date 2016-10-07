function D=roundDistanceMatrix2(D,tape,level)
% D=roundDistanceMatrix(D,tape,level)
% 
% tape - vector of margins
% level - vector of grid steps 
%
% for each distance d in matrix D:
%
%  if d > tape  it's coarsed infinity
%  if tape(i-1) < d < tape(i) it's coarsed to grid level(i)
% 


T=length(tape);

% assert(T==length(level));

D(D>tape(T)) = inf;

tape=[0 tape];
level=[0 level];

for t=T+1:-1:2
    I = find( D<=tape(t) & D>tape(t-1)  );
    D(I) = round( D(I)/level(t))*level(t);
end
