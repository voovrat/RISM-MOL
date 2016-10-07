function r_cs = rMHNC_closure(r,rGamma_in, EXPONENTA, u_short)
% !! size(r) = NxNumAtom*2

r_cs =  r.*exp(rGamma_in./r).*EXPONENTA - r - rGamma_in;

gamma = rGamma_in./r;
%c2 = gamma.*log(EXPONENTA) - gamma - 1;

g0 = (rGamma_in+r_cs)./r;
%g2 = gamma+c2;
   
I = find( (-r.*u_short + rGamma_in > 0 ) & ( EXPONENTA~=0 ) );

if ~isempty(I)
    r_cs(I) = -r(I).*u_short(I);
end
